#!/usr/bin/env python3
# -*- coding: utf-8
import os
import sys
import Bio.SeqIO
import torch
import numpy as np
from collections import defaultdict
from knn import BlastKnn, NetworkKnn
from utils import GeneOntology, get_term_scores, NS_ID, get_data, log
from network import read_network


pj_dir = os.path.dirname(os.path.abspath(__file__))
ontology = GeneOntology(os.path.join(pj_dir, 'data/go-basic.obo'))
train_seqs, train_pro_anno = get_data(os.path.join(pj_dir, 'data/train/train.seq'),
                                      os.path.join(pj_dir, 'data/train/train.txt'), ontology)


def split_func_scores(func_scores):
    res = {ns: defaultdict(dict) for ns in NS_ID}
    for seq in func_scores:
        for func in func_scores[seq]:
            res[ontology[func].namespace][seq][func] = func_scores[seq][func]
    return res

def do_naive(seqs):
    naive_res = os.path.join(pj_dir, 'model/naive/naive.txt')
    term_scores = defaultdict(dict)
    with open(naive_res) as fp:
        for line in fp:
            go_term, score = line.split()
            for seq in seqs:
                term_scores[seq.id][go_term] = float(score)
    return split_func_scores(term_scores)


def do_blast_knn(seqs, path):
    blast_db = os.path.join(pj_dir, 'model/blastknn/blastdb')
    blast_knn = BlastKnn(n_jobs=40)
    blast_knn.fit(train_seqs, train_pro_anno, blast_db=blast_db)
    return split_func_scores(blast_knn.predict(seqs, ontology=ontology, tmp_path=path,
                                               out_file=os.path.join(path, 'blast.out')))


def do_net_knn(seqs, seq_fasta, path):
    net_knn_db = os.path.join(pj_dir, 'model/netknn/networkdb/networkdb')
    network = read_network(os.path.join(pj_dir, 'model/netknn/network_cafa5.v11.5.txt'), 1.0)
    net_knn = NetworkKnn(network, n_jobs=40)
    net_knn.fit(train_seqs, train_pro_anno)
    sim_seqs_para = {'seqs_file': seq_fasta,
                     'db': net_knn_db,
                     'out_file': os.path.join(path, 'n_k_blast.out'),
                     'num_threads': 64}
    return split_func_scores(net_knn.predict(seqs, ontology=ontology, sim_seqs_para=sim_seqs_para))

# Requires modification
# def get_consensus_score(diff_func_scores, type_list, label_list=None, alpha=0.9):
#     res = defaultdict(dict)
#     for ns in NS_ID:
#         for seq in type_list[ns]:
#             term_list = label_list[ns] if label_list is not None \
#                 else reduce(lambda x, y: x | set(y[ns].get(seq, {})), diff_func_scores, set())
#             for go_term in term_list:
#                 res[seq][go_term] = 1.0
#                 for term_scores in diff_func_scores:
#                     res[seq][go_term] *= (1.0 - alpha * term_scores[ns].get(seq, {}).get(go_term, 0.0))
#                 res[seq][go_term] = 1.0 - res[seq][go_term]
#     return res
    
def main(argv):
    path = argv[0]
    log(path)
    seq_fasta = os.path.join(path, 'seqs.fasta')
    seqs = list(Bio.SeqIO.parse(seq_fasta, 'fasta'))
    
    naive = do_naive(seqs)
    log('naive')

    blast = do_blast_knn(seqs, path)
    log('blast_knn')

    net = do_net_knn(seqs, seq_fasta, path)
    log('net_knn')

    component_dict = {'naive': naive, 'blast': blast, 'net': net}
    for cop in component_dict:
        cop_res = component_dict[cop]
        with open(os.path.join(path, '%s_res.txt' % cop), 'w') as fcr:  # output component results
            for ns in cop_res:
                for seq in seqs:
                    scores = ontology.transfer_scores(cop_res[ns][seq.id])
                    for func, score in sorted(scores.items(), key=lambda x: x[1], reverse=True)[:100]:
                        print(seq.id, func, round(score, 3), sep='\t', file=fcr)


if __name__ == '__main__':
    log('begin')
    main(sys.argv[1:])
    log('finished')

# Usage: python run.py <path_to_query_folder>