#!/usr/bin/env python3
# -*- coding: utf-8
import sys
import os
import hashlib
import configparser
from collections import defaultdict
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast import NCBIXML
from functools import reduce
from scipy.sparse import csr_matrix
import time

from knn import Knn
from utils import get_data, output_res, GeneOntology, get_label_list, write_fasta

__all__ = ['BlastKnn', 'psiblast']


def psiblast(seqs, blast_db=None, evalue=0.001, out_file=None, num_iterations=1,
             tmp_path=None, num_threads=40, **kwargs):
    if not os.path.exists(out_file):
        if isinstance(seqs, str):
            seqs_file = seqs
        else:
            seqs_file = os.path.join(tmp_path, hashlib.md5(str(time.time()).encode()).hexdigest() + '.fasta')
            write_fasta(seqs, seqs_file)
        NcbipsiblastCommandline(query=seqs_file, db=blast_db, evalue=evalue, outfmt=5, out=out_file,
                                num_iterations=num_iterations, num_threads=num_threads, **kwargs)()
    with open(out_file) as fp:
        dis = []
        for rec in NCBIXML.parse(fp):
            query, sim = rec.query, []
            for alignment in rec.alignments:
                score = float('-inf')
                for hsp in alignment.hsps:
                    score = max(score, hsp.bits)
                    # score = max(score, -math.log10(hsp.expect) if hsp.expect > 0 else 250)
                pid = alignment.hit_def.split()[0]
                """If we want to run blast over training set, we should ignore itself."""
                if pid != query:
                    sim.append((pid, score))
            dis.append(sim)
    return dis


class BlastKnn(Knn):
    """

    """

    def __init__(self, **kwargs):
        self.blast_db = None
        super(BlastKnn, self).__init__(**kwargs)

    def fit(self, seqs, pro_anno, blast_db=None, **kwargs):
        self.blast_db = blast_db
        super(BlastKnn, self).fit(seqs, pro_anno, **kwargs)

    def similarity(self, seqs, **kwargs):
        """
        :param seqs: protein seqs
        :param kwargs: psiblast parameters
        :return: psiblast score
        """
        # psiblast(seqs, **kwargs)
        # if self.n_jobs > 1:
        #     pool = Pool(self.n_jobs)
        #     res, size = [],  len(seqs) // self.n_jobs + (len(seqs) % self.n_jobs > 0)
        #     for i in range(min(self.n_jobs, len(seqs))):
        #         res.append(pool.apply_async(psiblast, (seqs[i*size:min((i+1)*size, len(seqs))],), kwargs))
        #     pool.close()
        #     pool.join()
        #     return list(itertools.chain(*[r.get() for r in res]))
        dist = psiblast(seqs, blast_db=self.blast_db, num_threads=self.n_jobs, **kwargs)
        row, col, data = [], [], []
        for i in range(len(dist)):
            for pid, score in dist[i]:
                row.append(i)
                col.append(self.train_pid[pid])
                data.append(score)
        return csr_matrix((data, (row, col)), shape=(len(dist), self.trainy.shape[0]))


def main(argv):
    if len(argv) == 0:
        print('please input one or more configure files')
        print('train.seqs           ', 'the sequences FASTA file of training set')
        print('train.annotations    ', 'the annotations file of training set')
        print('benchmark.seqs       ', 'the sequences FASTA file of benchmarks')
        print('output.results       ', 'the output file path')

    conf = configparser.ConfigParser()
    print(conf.read(argv))

    ontology = GeneOntology(conf.get('train', 'ontology'))
    blast_knn = BlastKnn(n_jobs=conf.getint('model', 'n_jobs'))

    train_seqs, train_pro_anno = get_data(conf.get('train', 'seqs'), conf.get('train', 'annotations'), ontology)
    label_list = get_label_list(conf.get('label', 'label_list', fallback=None), ontology)
    if label_list is not None:
        label_list = reduce(lambda x, y: x | y, label_list.values())
    if conf.has_section('target'):
        species = conf.get('target', 'species')
        print(species)
        test_seqs, _ = get_data(conf.get('target', 'seqs') + species + '.fasta')
        blast_out, res_out = species + '.xml', species + '.res'
    else:
        test_seqs, _ = get_data(conf.get('fasta', 'all', fallback=conf.get('train', 'seqs')), None)
        blast_out = res_out = ''
    blast_knn.fit(train_seqs, train_pro_anno, label_list=label_list, blast_db=conf.get('model', 'blastdb_path'))
    res = blast_knn.predict(test_seqs, num_iterations=conf.getint('model', 'num_iterations', fallback=1),
                            ontology=ontology if conf.getboolean('model', 'normalization', fallback=True) else None,
                            out_file=conf.get('model', 'out') + blast_out,
                            tmp_path=conf.get('model', 'tmp_path'),
                            normalized=conf.getboolean('model', 'normalized', fallback=True))
    output_res(conf.get('output', 'results') + res_out, res, keywords='Blast-Knn')


if __name__ == '__main__':
    main(sys.argv[1:])
