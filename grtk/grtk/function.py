import os, gzip, urllib, itertools

from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
import mysql.connector
import numpy
import pandas

import grtk.bbi

def queryUCSC(assembly, query):
    db = mysql.connector.connect(user="genome", password="genome",
                                 host="wrendb", database=assembly)
    df = pandas.io.sql.read_frame(query, db)
    db.close()
    return df

def open_url(url):
    """Download the file at the given URL if it hasn't been previously downloaded,
    and return a handle."""
    name = url.split("/")[-1]
    path = "/tmp/" + name
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path) 
    if path.endswith(".gz"):
        return gzip.open(path)
    else:
        return open(path)

def gene2go(taxon_id):
    taxon_id = int(taxon_id)
    mapping = {}
    with open_url("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz") as h:
        h.readline()
        for line in h:
            taxon, gene_id, go_id = line.split()[:3]
            if taxon_id != int(taxon):
                continue
            gene_id = int(gene_id)
            go_id = int(go_id[3:])
            if not gene_id in mapping:
                mapping[gene_id] = set()
            mapping[gene_id].add(go_id)
    return mapping
    

def read_matrix(bw_set, regions):
    X = {}
    for id, group in regions.groupby("name").groups.items():
        x = [bw_set.mean(chrom,start,end) 
             for chrom, start, end 
             in genes.ix[group,["chrom","txStart","txEnd"]].values]
        x.sort(key=lambda v: v.mean())
        X[id] = x[-1]
    return pandas.DataFrame(X)

# genes = queryUCSC("hg19", """
#      SELECT chrom, txStart, txEnd, 
#         knownToLocusLink.value as name, 0 as score, strand, 
#         cdsStart, cdsEnd, exonCount, exonStarts, exonEnds
#      FROM knownGene
#      INNER JOIN knownToLocusLink
#      ON knownToLocusLink.name=knownGene.name""")

# lincRNAs = queryUCSC("hg19", """SELECT * from lincRNAsTranscripts""")

# bws = grtk.bbi.BigWigSet("/home/gilesc/data/continuous-tracks/hg19/RNAseq/")

# go_human = gene2go(9606)
# X = read_matrix(bws, genes)


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
        

def train(df, annotations):
    X = df.as_matrix().T
    Y = [list(annotations.get(int(gene),set())) for gene in df.columns]
    #inner = SVC(kernel="linear", probability=True)
    #inner = LogisticRegression()
    inner = SimplePrior()
    clf = OneVsRestClassifier(inner)
    clf.fit(X, Y)
    return clf

def evaluate(X, model, annotations):
    loss = []
    for gene_id in X.columns:
        M = annotations.get(int(gene_id), set())
        y = numpy.array([1 if label in M else 0 for label in model.classes_])
        y_hat = [[p, 1-p] for p in model.predict_proba(X[gene_id])[0]]
        loss.append(log_loss(y, y_hat))
    return numpy.array(loss)
