#!/usr/bin/env python3
# -*- coding: utf-8
"""
Created on 2016/9/5
@author yrh

"""

import sys
import configparser
from collections import defaultdict
sys.path.append("/home/zhushanfeng/storage/liuhc/baseline/PFP")
from utils import get_data, output_res, GeneOntology, NS_ID




def main(argv):
    if len(argv) == 0:
        print('please input one or more configure files')
        print('train.seqs           ', 'the sequences FASTA file of training set')
        print('train.annotations    ', 'the annotations file of training set')
        print('benchmark.seqs       ', 'the sequences FASTA file of benchmarks')
        print('output.results       ', 'the output file path')

    conf = configparser.ConfigParser()
    print(conf.read(argv))
    ontology = GeneOntology(conf.get('train', 'ontology', fallback=None))
    with open(conf.get('output', 'results'), 'w') as fout:
        for ns in NS_ID:
            _, pro_anno = get_data(conf.get('train', 'seqs'), conf.get('train', ns), ontology)
            total = len(pro_anno)
            anno_count = defaultdict(int)
            for go_terms in pro_anno.values():
                for go_term in go_terms:
                    anno_count[go_term] = anno_count[go_term] + 1 / total
            pred = sorted(anno_count.items(), key=lambda x: x[1], reverse=True)
            for item in pred:
                print(item[0], item[1], file=fout)


if __name__ == '__main__':
    main(sys.argv[1:])
