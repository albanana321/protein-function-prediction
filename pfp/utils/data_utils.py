#!/usr/bin/env python3
# -*- coding: utf-8
import os
import sys
import datetime
import configparser
from Bio import SeqIO
from collections import defaultdict
import json

__all__ = ['Annotation',
           'parse_fasta',
           'write_fasta',
           'get_data',
           'get_target_data',
           'generate_target_data',
           'output_res',
           'split_annotations_by_namespace',
           'get_protein_species',
           'get_type_list',
           'get_blast_identity',
           'get_sample',
           'get_now',
           'NS_ID',
           'NS_ID_FC2S',
           'get_term_scores',
           'log']


NS_ID = ['mf', 'bp', 'cc']
NS_ID_FC2S = {'molecular_function': 'mf', 'biological_process': 'bp', 'cellular_component': 'cc',
              'F': 'mf', 'P': 'bp', 'C': 'cc'}
result_key_words = ('AUTHOR', 'MODEL', 'KEYWORDS', 'ACCURACY', 'END')
root_term = ('GO:0008150', 'GO:0003674', 'GO:0005575')


class Annotation(defaultdict):
    """
    Protein Annotations.
    """

    def __init__(self):
        super(Annotation, self).__init__(set)

    def combine(self, *args):
        """
        :param args: one or more ProteinAnno
        :return: combined ProteinAnno
        """
        for pro_anno in args:
            for name in pro_anno:
                self[name] = self[name] | pro_anno[name]
        return self

    @staticmethod
    def load(anno_files=None, ontology=None):
        """
        :param anno_files: the BP, MF, CC functions files, the file have two columns: ID  GO
        :param ontology: if it is not None, the GO will propagate with it
        :return: instance of ProteinAnno
        """
        if anno_files is None:
            return None
        if isinstance(anno_files, str):
            anno_files = [anno_files]
        pro_anno = Annotation()
        for funcs_file in anno_files:
            try:
                with open(funcs_file) as fp:
                    for line in fp:
                        pid, go_term = line.split()[:2]
                        pro_anno[pid].add(go_term)
            except Exception:
                pass
        if ontology is not None:
            for pid in pro_anno:
                pro_anno[pid] = ontology.transfer(pro_anno[pid])
        else:
            print('Warning! Don\'t transfer the annotations!')
        return pro_anno


def get_term_scores(res_file, ontology):
    if res_file is None:
        return None
    term_scores = defaultdict(dict)
    for ns in NS_ID:
        term_scores[ns] = defaultdict(dict)
    with open(res_file) as fp:
        for line in fp:
            try:
                if line.startswith(result_key_words):
                    continue
                line_list = line.split()
                pid, go_term, score, *_ = line_list
                if go_term in ontology:
                    term_scores[ontology[go_term].ns][pid][go_term] = float(score)
                    term_scores[pid][go_term] = float(score)
            except Exception:
                pass
        return term_scores


def parse_fasta(fasta_file):
    return SeqIO.parse(fasta_file, 'fasta')


def write_fasta(seqs, fasta_file):
    SeqIO.write(seqs, fasta_file, 'fasta')


def get_data(train_seq_file, train_anno_file=None, ontology=None):
    return list(parse_fasta(train_seq_file)), Annotation.load(train_anno_file, ontology)


def get_target_data(targets_file, targets_dir):
    if not os.path.exists(targets_file):
        seqs = generate_target_data(targets_dir)
        write_fasta(seqs, targets_file)
    return parse_fasta(targets_file)


def generate_target_data(targets_dir):
    seqs = set()
    for dirpath, dirnames, filenames in os.walk(targets_dir):
        for file in filenames:
            if file.endswith('.tfa'):
                seqs |= set(parse_fasta(os.path.join(dirpath, file)))
        for dname in dirnames:
            seqs |= generate_target_data(os.path.join(dirpath, dname))
    return seqs


def output_res(res_file, pred, author='FDUPFP', model_id=1, keywords='machine learning'):
    if res_file is None: return
    with open(res_file, 'w') as fp:
        # print('AUTHOR', author, file=fp)
        # print('MODEL', model_id, file=fp)
        # print('KEYWORDS', keywords + '.', file=fp)
        for seq in pred:
            for func, score in sorted(pred[seq].items(), key=lambda x: x[1], reverse=True)[:200]:
                print(seq, func, round(score, 3), sep='\t', file=fp)
        # print('END', file=fp)


def split_annotations_by_namespace(data_file, bp_file, cc_file, mf_file):
    annotations = {}
    with open(data_file) as fp:
        for line in fp:
            pid, acc, np = line.strip().split()[:3]
            if pid not in annotations:
                annotations[pid] = {'bp': set(), 'cc': set(), 'mf': set()}
            annotations[pid][np].add(acc)
    with open(bp_file, 'w') as fbp, open(cc_file, 'w') as fcc, open(mf_file, 'w') as fmf:
        fp = {'bp': fbp, 'cc': fcc, 'mf': fmf}
        for pid in annotations:
            for np in annotations[pid]:
                if annotations[pid][np] == {'GO:0005515'}:
                    continue
                for acc in annotations[pid][np]:
                    print(pid, acc, file=fp[np])


def get_protein_species(protein_species_file):
    protein_species = {}
    with open(protein_species_file) as fp:
        for line in fp:
            pid, sp = line.split()[:2]
            protein_species[pid] = sp
    return protein_species


def get_type_list(type_file):
    if type_file is None:
        return None
    type_list = set()
    with open(type_file) as fp:
        for line in fp:
            type_list.add(line.strip())
    return type_list


def get_blast_identity(file, binary=None):
    if file is None:
        return None
    blast_id = {}
    with open(file) as fp:
        for line in fp:
            pid, id = line.split()
            blast_id[pid] = float(id) if binary is None else int(float(id) >= binary)
    return blast_id

def get_sample(sample_file, times=10000):
    with open(sample_file) as fp:
        samples = [line.strip() for line in fp]
    length = len(samples) // times
    return [samples[i*length: (i+1)*length] for i in range(times)]

def get_now():
    return datetime.datetime.now().strftime('%y-%m-%d %H:%M:%S')

def log(*args):
    print(f'[{datetime.datetime.now()}]', *args)

def main(argv):
    conf = configparser.ConfigParser()
    print(conf.read(argv))
    split_annotations_by_namespace(conf.get('annotations', 'data'),
                                   conf.get('annotations', 'bp'),
                                   conf.get('annotations', 'cc'),
                                   conf.get('annotations', 'mf'))


if __name__ == '__main__':
    main(sys.argv[1:])

