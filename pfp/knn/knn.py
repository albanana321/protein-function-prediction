#!/usr/bin/env python3
# -*- coding: utf-8
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from scipy.sparse import csc_matrix
from utils import save, load, log


__all__ = ['Knn']

class BaseModel(object, metaclass=ABCMeta):
    """
    Based Classifier.
    """
    def __init__(self, model_path=None, n_jobs=1):
        """
        :param model_path: the path to save the model.
        :param n_jobs: the parallel jobs.
        """
        self.model_path = model_path
        self.n_jobs = int(n_jobs)

    @abstractmethod
    def fit(self, X, y, **kwargs):
        pass

    @abstractmethod
    def predict(self, X, **kwargs):
        pass

    @staticmethod
    def save(save_path, content, key='model'):
        save(save_path, content, key)

    @staticmethod
    def load(load_path, key='model'):
        return load(load_path, key)
    
class Knn(BaseModel, metaclass=ABCMeta):
    """

    """

    def __init__(self, **kwargs):
        self.go2id = self.id2go = self.go_mat = self.seqs = self.trainy = self.train_pid = None
        super(Knn, self).__init__(**kwargs)

    def fit(self, seqs, pro_anno, label_list=None, **kwargs):
        self.go2id, self.id2go, self.train_pid, row, col, data = {}, [], {}, [], [], []
        for idx, seq in enumerate(seqs):
            self.train_pid[seq.id] = idx
            for go_term in pro_anno[seq.id]:
                if label_list is not None and go_term not in label_list:
                    continue
                if go_term not in self.go2id:
                    self.go2id[go_term] = len(self.go2id)
                    self.id2go.append(go_term)
                row.append(idx)
                col.append(self.go2id[go_term])
                data.append(1)
        self.trainy = csc_matrix((data, (row, col)), shape=(len(seqs), len(self.id2go)))
        log('fit finish')

    def predict(self, seqs, ontology=None, normalized=True, **kwargs):
        score_mat = self.similarity(seqs, **kwargs).dot(self.trainy)
        term_scores = defaultdict(dict)
        for idx, seq in enumerate(seqs):
            max_ns, scores = defaultdict(float), {}
            for col, score in zip(score_mat[idx].indices, score_mat[idx].data):
                go_term = self.id2go[col]
                scores[go_term] = score
                if ontology is not None and normalized:
                    ns = ontology[go_term].ns
                    max_ns[ns] = max(max_ns[ns], score)
            if ontology is None or not normalized:
                max_ns = defaultdict(lambda: 1.0)
            for go_term in scores:
                term_scores[seq.id][go_term] = scores[go_term] / max_ns[ontology[go_term].ns]
        return term_scores

    @abstractmethod
    def similarity(self, seqs, **kwargs):
        pass
