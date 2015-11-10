#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tempfile
import os
import itertools
from collections import defaultdict
from subprocess import Popen, PIPE, check_output
from os.path import realpath, dirname

import numpy
from pandas import DataFrame, Series, SparseDataFrame
import pandas
from sklearn.preprocessing import Imputer
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

import grtk.db

db = grtk.db.connect()
c = db.cursor()
DataFrame.show = lambda self: self.head(n=3).T.head(n=3).T
SparseDataFrame.show = DataFrame.show

def which(exe):
    return realpath(check_output(["which", exe]).strip())

def multimap(pairs):
    m = {}
    for k,v in pairs:
        if not k in m:
            m[k] = set()
        m[k].add(v)
    for k,vs in m.items():
        m[k] = tuple(vs)
    return m

def impute(X):
    """
    Impute the value of missing genes to be the mean expression of that gene
    across all experiments.
    
    Also drops experiments and genes with no data.
    """
    if not pandas.isnull(X).sum().sum():
        return X
    X = X.dropna(axis=0, how="all").dropna(axis=1, how="all")
    model = Imputer(axis=1, strategy="mean")
    return DataFrame(model.fit_transform(X),
                     index=X.index,
                     columns=X.columns)

URSA_MEAN = 7.32
URSA_STDEV = 2.63

def symbol_map():
    m = {}
    c.execute("""
    SELECT id, symbol 
    FROM gene 
    WHERE taxon_id=9606 AND symbol IS NOT NULL""")
    return dict(c)

def collapse_to_symbols(X, axis=1, min_pct=0.2):
    # FIXME: collapse by max mean instead of max
    thresh = int(X.shape[0] * min_pct)
    return X.groupby(symbol_map(), axis=1).max().dropna(thresh=thresh, axis=1)

def fetch_genes(taxon_id):
    c.execute("""
    SELECT id, symbol, name 
    FROM gene 
    WHERE taxon_id=%s 
    ORDER BY id""", (taxon_id,))
    return DataFrame.from_items([(row[0], row) for row in c], 
                                columns=["id", "symbol", "name"],
                                orient="index")

def platform_expression(accession, require_age=True, require_gender=False, limit=None):
    genes = fetch_genes(9606)
    query = """
    SELECT sample.id, sample.age, sample.gender, expression.data 
    FROM expression 
    INNER JOIN sample 
    ON expression.sample_id=sample.id 
    INNER JOIN platform
    ON sample.platform_id=platform.id
    WHERE platform.accession=%s"""
    if require_age:
        query += "\nAND sample.age IS NOT NULL"
    if require_gender:
        query += "\nAND sample.gender IS NOT NULL"
    if limit:
        query += "\tLIMIT " + str(limit)
    c.execute(query, (accession,))
    samples, age, gender, expression = zip(*c)
    X = DataFrame.from_records(list(expression), 
                               index=samples, columns=genes.index)
    X.index.name = "Sample ID"
    X.columns.name = "Gene ID"
    P = DataFrame({"age": age, "gender": gender}, index=samples)
    P.index.name = "Sample"
    return X, P

def sample_expression(accession):
    c.execute("""
    SELECT data FROM expression
    INNER JOIN sample
    ON sample.id=expression.sample_id
    WHERE sample.accession=%s""", (accession,))
    try:
        return Series(c.__next__()[0], 
                      index=fetch_genes(9606).index)
    except StopIteration:
        return None

def infer_tissue(X):
    """
    Infer tissue for sample(s). The columns should be gene symbols.
    """
    if X.columns.dtype != 'O':
        raise Exception("The expression matrix columns are not object types, implying they do not contain gene symbols. Gene symbols are required for tissue imputation.")

    URSA_BIN = which("ursa")
    Xt = impute(X).T
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "24"
    with tempfile.NamedTemporaryFile() as h:
        X.T.dropna().to_csv(h.name, sep="\t")
        p = Popen([URSA_BIN, "-i", h.name, "-n", str(X.shape[0]), "-r", "1"],
                  cwd=dirname(URSA_BIN), stdout=PIPE, env=env)
        T = pandas.io.parsers.read_csv(p.stdout, skiprows=(6*X.shape[0] + 1),
                                       sep="\t", index_col=0, header=False,
                                       names=X.index)
    c.execute("""
    SELECT translate(term.name, ' ', '_'), term.id
    FROM term
    INNER JOIN ontology
    ON ontology.id=term.ontology_id
    WHERE ontology.namespace='BTO'""")
    m = dict([(k.replace("/", "--"),v) for k,v in c])
    T.index = [m[k] for k in T.index]
    T.index.name = "Term ID"
    T.columns.name = "Sample ID"
    return T

def insert_sample_terms(T, min_p=0.5):
    """
    Insert sample-term associations.
    """
    assert(T.index.name=="Term ID")
    assert(T.columns.name=="Sample ID")
    #for sample_id in T.index:
    #    c.execute("""
    #    DELETE FROM sample_term 
    #    INNER JOIN sample 
    #    ON sample.id=sample_term.sample_id
    #    WHERE sample.id=%s""", (sample_id,))
    c.executemany("""
    INSERT INTO sample_term (sample_id, term_id, probability)
    VALUES (%s,%s,%s)""", 
                  ((int(r[1]), int(r[2]), float(r[3])) for r in
                   T[T>min_p].unstack().reset_index().dropna().to_records()))
    db.commit()


def coo_to_df(triplets):
    """
    Create a SparseDataFrame from a sequence of (row,col,value) triplets.
    """
    data = defaultdict(dict)
    for row, col, val in triplets:
        data[col][row] = val
    return SparseDataFrame(data)

def tissue_expression_training_set(taxon_id=9606, limit=200):
    c.execute("""
    SELECT sample_term.sample_id, expression.data, 
        sample_term.term_id, sample_term.probability
    FROM sample_term
    INNER JOIN term
    ON term.id=sample_term.term_id
    INNER JOIN ontology
    ON ontology.id=term.ontology_id
    INNER JOIN sample
    ON sample.id=sample_term.sample_id
    INNER JOIN expression
    ON expression.sample_id=sample.id
    INNER JOIN platform
    ON sample.platform_id=platform.id
    INNER JOIN taxon
    ON platform.taxon_id=taxon.id
    WHERE ontology.namespace='BTO'
    AND sample_term.probability=1
    AND taxon.id=%s
    ORDER BY random()
    LIMIT %s""", (taxon_id, limit))
    samples, data, tissues, values = zip(*c)
    T = coo_to_df(zip(samples, tissues, values))
    T.index.name = "Sample ID"
    T.columns.name = "Term ID"
    c.execute("""SELECT id FROM gene WHERE gene.taxon_id=%s ORDER BY id""", 
              (taxon_id,))
    X = DataFrame.from_records(list(data),
                               index=samples, columns=[e[0] for e in c])
    return X,T

def ontology_graph():
    import networkx as nx
    c.execute("""
    SELECT agent, target 
    FROM term_relation
    WHERE probability=1""")
    g = nx.DiGraph()
    g.add_edges_from(c)
    return g


# Multilabel learning

class BRRDT(object):
    """
    Implements 
    https://www.siam.org/proceedings/datamining/2010/dm10_068_zhangx.pdf
    """
    def __init__(self, n=10):
        self._trees = []

    def fit(self, X, Y):
        pass

from sklearn.base import BaseEstimator, ClassifierMixin
class SimplePrior(BaseEstimator, ClassifierMixin):
    "A very simple classifier that always predicts the prior probability."
    def __init__(self):
        pass
    
    def fit(self, X, Y):
        p = Y.sum() / Y.size
        self.p = numpy.array([p, 1-p])
        
    def predict_proba(self, X):
        return numpy.array([self.p,] * X.shape[0])
        

def train(X, Y):
    """
    Train a multilabel classifier.
    - X is a DataFrame of instances (rows) versus features (columns).
    - Y is a SparseDataFrame of instances (rows) versus term IDs (columns).
    """
    #X = impute(X)
    #common = list(set(X.index) & set(Y.index))
    #X = X.ix[common,:]
    #Y = Y.ix[common,:]
    Yt = Y.T
    Yl = []
    for sample in Yt.columns:
        Yl.append(tuple(Yt.ix[:,sample].to_dense().notnull().nonzero()[0]))
    #inner = SVC(kernel="linear", probability=True)
    inner = LogisticRegression()
    #inner = SimplePrior()
    clf = OneVsRestClassifier(inner)#, n_jobs=-1)
    clf.fit(X, Yl)
    return clf

def evaluate(X, model, annotations):
    loss = []
    for gene_id in X.columns:
        M = annotations.get(int(gene_id), set())
        y = numpy.array([1 if label in M else 0 for label in model.classes_])
        y_hat = [[p, 1-p] for p in model.predict_proba(X[gene_id])[0]]
        loss.append(log_loss(y, y_hat))
    return numpy.array(loss)

#C = X.corrwith(P)
#C = C.ix[numpy.invert(numpy.isnan(C))]
#C.sort()
#db.close()

class SparseMultilabelLR(object):
    def __init__(self):
        pass
