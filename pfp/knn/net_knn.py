#!/usr/bin/env python3
# -*- coding: utf-8
import os
import sys
import configparser
from functools import reduce
from scipy.sparse import csr_matrix
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Blast import NCBIXML

from knn import Knn
from utils import GeneOntology, get_data, output_res, get_label_list
from network import read_network


__all__ = ['NetworkKnn', 'get_sim_seqs']


def get_sim_seqs(seqs_file=None, db=None, evalue=0.01, out_file=None, num_iterations=1, **kwargs):
    if out_file is None:
        return None
    if not os.path.exists(out_file):
        NcbipsiblastCommandline(query=seqs_file, db=db, evalue=evalue, outfmt=5, out=out_file,
                                num_iterations=num_iterations, **kwargs)()
    sim_seqs = {}
    with open(out_file) as fp:
        for rec in NCBIXML.parse(fp):
            max_score = 0, float('-inf'), None
            for alignment in rec.alignments:
                score = float('-inf')
                for hsp in alignment.hsps:
                    score = max(score, hsp.bits)
                pid = alignment.hit_def.split()[0]
                # if pid == rec.query:
                #     continue
                max_score = max(max_score, (pid == rec.query, score, pid))
            sim_seqs[rec.query] = max_score[2], max_score[1]
    return sim_seqs


class NetworkKnn(Knn):
    """

    """

    def __init__(self, network, **kwargs):
        self.network, self.train_id = network, {}
        super(NetworkKnn, self).__init__(**kwargs)

    def fit(self, seqs, pro_anno, **kwargs):
        tot, net_seqs = 0, []
        for seq in seqs:
            if seq.id in self.network:
                self.train_id[seq.id] = tot
                net_seqs.append(seq)
                tot += 1
        super(NetworkKnn, self).fit(net_seqs, pro_anno, **kwargs)

    def similarity(self, seqs, n_neighbors=50, sim_seqs=None, **kwargs):
        if sim_seqs is None:
            sim_seqs = {}
        row, col, data = [], [], []
        for idx, seq in enumerate(seqs):
            if seq.id not in self.network and seq.id not in sim_seqs:
                continue
            sim_seq_id = sim_seqs[seq.id][0]
            #assert (seq.id in self.network and sim_seq_id == seq.id) or seq.id not in self.network
            dist = []
            for pid in self.network[sim_seq_id]:
                if pid in self.train_id:
                    dist.append((pid, self.network[sim_seq_id][pid]))
            dist = sorted(dist, key=lambda x: x[1], reverse=True)
            for pid, weight in dist[:n_neighbors]:
                row.append(idx)
                col.append(self.train_id[pid])
                data.append(weight)
        return csr_matrix((data, (row, col)), shape=(len(seqs), self.trainy.shape[0]))

    def predict(self, seqs, sim_seqs=None, sim_seqs_para=None, **kwargs):
        if sim_seqs is None:
            sim_seqs = get_sim_seqs(**sim_seqs_para)
        return super(NetworkKnn, self).predict(seqs, sim_seqs=sim_seqs, **kwargs)


def main(argv):
    conf = configparser.ConfigParser()
    print(conf.read(argv))

    print(conf.get('network', 'output'))
    ontology = GeneOntology(conf.get('train', 'ontology'))
    network = read_network(conf.get('network', 'output'), conf.getfloat('network', 'self2self', fallback=1.0))
    sim_seqs_para = {'seqs_file': conf.get('fasta', 'all'),
                     'db': conf.get('model', 'blastdb_path'),
                     'out_file': conf.get('model', 'out'),
                     'num_threads': conf.get('model', 'n_jobs', fallback=200)}

    train_seqs, train_pro_anno = get_data(conf.get('train', 'seqs'), conf.get('train', 'annotations'), ontology)
    label_list = get_label_list(conf.get('label', 'label_list', fallback=None), ontology)
    if label_list is not None:
        label_list = reduce(lambda x, y: x | y, label_list.values())

    model = NetworkKnn(network)
    model.fit(train_seqs, train_pro_anno, label_list=label_list)
    print('finish fit')
    test_seqs, _ = get_data(conf.get('fasta', 'all'))
    res = model.predict(test_seqs, ontology=ontology, sim_seqs_para=sim_seqs_para,
                        n_neighbors=conf.getint('model', 'n_neighbors', fallback=50),
                        normalized=conf.getboolean('model', 'normalized', fallback=True))
    output_res(conf.get('output', 'results', fallback=None), res, keywords='Network-Knn')


if __name__ == '__main__':
    main(sys.argv[1:])
