from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.svm import SVR
from numpy import array, fabs
import numpy as np

import imp
import gfa
imp.reload(gfa)

class SimpleRegressor(object):
    def __init__(self, method="mean"):
        assert(method in ("mean", "mode"))
        self._method = method

    def fit(self, X, y):
        if self._method == "mean":
            self._y = y.mean()
        elif self._method == "mode":
            self._y = np.argmax(np.bincount(np.round(y).astype(int)))
    
    def predict(self, X):
        return array([self._y for _ in range(X.shape[0])])

def impute_age():
    X, P = gfa.platform_expression("GPL96")
    model = impute.KNNImputer()
    Xi = model.fit_transform(X, axis=1)

    age = array(P["age"].tolist())
    Xm = Xi.as_matrix()
    ix = array((age >= 10) & (age <= 120)).nonzero()[0]
    np.random.shuffle(ix)
    Xm = Xm[ix,:]
    age = age[ix]

    n_train = 2000
    n_test = 500
    #clf = SVR(C=1e-5, epsilon=1)
    #clf = LinearRegression()
    clf = Ridge()
    #clf = SimpleRegressor()
    #clf = Lasso()
    clf.fit(Xm[:n_train,:], age[:n_train])
    y = age[n_train:(n_train+n_test)]
    y_hat = clf.predict(Xm[n_train:(n_train+n_test)])
    dy = y - y_hat

    bias_tr = y_hat.mean() - age.mean()
    print("\nBias (vs train):\t\t", bias_tr)
    print("Bias (vs test):\t\t\t", dy.mean())
    print("Mean error:\t\t\t", fabs(dy).mean())
    print("Mean error (bias corrected):\t", fabs(dy - bias_tr).mean())
    print("MSE:\t\t\t\t", np.power(dy,2).mean())
    #clf.fit(Xi.as_matrix(), age)

from sklearn.datasets import make_multilabel_classification
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelBinarizer, Imputer
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC, SVC
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest
from pandas import DataFrame, SparseDataFrame

from grtk.impute import KNNImputer
 
"""
g = gfa.ontology_graph()
X,T = gfa.tissue_expression_training_set(limit=10000)
X = X.dropna(axis=1, how="all")
#Xi = KNNImputer().fit_transform(X)
Xi = DataFrame(Imputer().fit_transform(X.as_matrix()),
               index=X.index, columns=X.columns)
Xs = Xi.as_matrix()

Tt = T.T
data = {}
for i in Tt.columns:
    data[i] = {}
    for j in Tt.index[Tt.ix[:,i].to_dense().notnull()]:
        data[i][j] = 1
        if j in g:
            for ancestor in g.successors(j):
                data[i][ancestor] = 1

Ts = SparseDataFrame(data).T.to_sparse(fill_value=0)
Ts = Ts.fillna(-1).as_matrix()

Xs = Xi.as_matrix()

base_clf = Pipeline([('feature_selection', SelectKBest(k=100)),
                     ('svm', SVC(probability=True))])
clf = OneVsRestClassifier(base_clf, n_jobs=-1)
N_train = 5000
ix = np.invert((Ts[:N_train,:] == -1).all(axis=0))
clf.fit(Xs[:N_train,:], Ts[:N_train,ix])

N_test = 1000
print(roc_auc_score(Ts[N_train:(N_train+N_test),ix].flat, 
                    clf.predict_proba(Xs[N_train:(N_train+N_test),:]).flat))
"""

import gfa

#X, P = gfa.platform_expression("GPL96")
#X = X.dropna(axis=1, how="all")
#Xi = DataFrame(Imputer().fit_transform(X.as_matrix()),
#               index=X.index, columns=X.columns)

#C = Xi.corrwith(P["age"])
#C.sort()

#ix = (P["age"] > 10) & (P["age"] < 100)
#C2 = Xi.ix[ix,:].corrwith(P.ix[ix,"age"])
#C2.sort()
