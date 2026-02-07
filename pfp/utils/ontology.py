#!/usr/bin/env python3
# -*- coding: utf-8

from collections import defaultdict
import sys
from goatools.obo_parser import GODag
from utils import NS_ID_FC2S


__all__ = ['GeneOntology']


class GOTerm(object):
    """

    """

    def __init__(self, go_term):
        """
        :param go_term: instance of the GOTerm
        :param godag: instance of the GODag
        :return:
        """
        self.id = go_term.id
        self.parents = {p.id for p in go_term.parents}
        if hasattr(go_term, 'relationship'):
            for parent in go_term.relationship.get('part_of', set()):
                if parent.namespace == go_term.namespace:
                    self.parents.add(parent.id)
        self.name = go_term.name
        self.namespace = self.ns = NS_ID_FC2S[go_term.namespace]
        self.children = set()
        self.depth = 0


class GeneOntology(dict):
    """

    """

    def __init__(self, obo_file_path):
        """
        :param obo_file_path: the ontology obo file
        :return:
        """
        super(GeneOntology, self).__init__()
        go_dag = GODag(obo_file_path, 'relationship')
        for go_id, go_term in go_dag.items():
            self[go_id] = GOTerm(go_term)
        self.get_children()
        self.root_term = ['GO:0008150', 'GO:0003674', 'GO:0005575']
        self.get_depth()

    def transfer(self, go_list):
        """
        :param go_list: the go terms which should be transferred
        :return:
        """
        go_list = list(filter(lambda x: x in self, go_list))
        ancestors, now = set(go_list), set(go_list)
        while len(now) > 0:
            next = set()
            for go_term in now:
                if go_term in self:
                    next |= self[go_term].parents - ancestors
            now = next
            ancestors |= now
        return ancestors

    def get_children(self):
        for go_id in self:
            for parent in self[go_id].parents:
                self[parent].children.add(go_id)

    def get_depth(self):
        for root in self.root_term:
            self[root].depth = 1
        now = set(self.root_term)
        while len(now) > 0:
            next = set()
            for go_term in now:
                for child in self[go_term].children:
                    if self[child].depth == 0:
                        next.add(child)
                        self[child].depth = self[go_term].depth + 1
            now = next

    def transfer_scores(self, term_scores):
        scores = defaultdict(float)
        for go_term in sorted(self.transfer(term_scores.keys()), key=lambda x: self[x].depth, reverse=True):
            scores[go_term] = max(scores[go_term], term_scores.get(go_term, 0))
            for parent_id in self[go_term].parents:
                scores[parent_id] = max(scores[parent_id], scores[go_term])
        return scores


if __name__ == '__main__':
    ontology = GeneOntology('go-basic.obo')
    print(ontology['GO:0090645'].id)
