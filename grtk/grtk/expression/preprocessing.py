import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
from ctypes import c_bool, c_double
import numpy as np
import pandas as pd

def standardize(X):
    """
    Standardize each row in X to mean = 0 and SD = 1.
    """
    X_m = np.ma.masked_invalid(X)
    return ((X.T - X_m.mean(axis=1)) / X_m.std(axis=1)).T.data

mask = None
X_s = None
X = None
k = None

def knn_init(k_, mask_, X_, X_s_):
    global k, mask, X_s, X
    mask = from_shared(mask_)
    X_s = from_shared(X_s_)
    X = from_shared(X_)
    k = k_

def knn_work(i):
    print(i)
    dx = X_s.dot(X_s[i,:]) / ((~ mask) & (~ mask[i,:])).sum(axis=1)
    ix = (-dx).argsort()
    for j in np.isnan(X[i,:]).nonzero()[0]:
        v = X[ix,j]
        v = v[np.invert(np.isnan(v))]
        X[i,j] = v[:k].mean()

def ctype_to_dtype(ctype):
    if ctype == c_double:
        return np.float64
    elif ctype == c_bool:
        return np.bool
    else:
        raise Exception

def to_shared(arr, type=c_double):
    shared = RawArray(type, arr.flat)
    return (shared, ctype_to_dtype(type), arr.shape)

def from_shared(args):
    arr, dtype, shape = args
    return np.frombuffer(arr, dtype=dtype).reshape(shape)

class KNNImputer(object):
    def __init__(self, k=50):
        self._k = k

    def fit_transform(self, X, axis=0):
        assert(axis in (0,1))
        if isinstance(X, pd.DataFrame):
            X = X.dropna(axis=0, how="all").dropna(axis=1, thresh=self._k)
            return pd.DataFrame(
                self.fit_transform(X.as_matrix(), axis=axis),
                index=X.index,
                columns=X.columns)
        if axis==0:
            return self.fit_transform(X.T, axis=1).T

        X_s = standardize(X)
        mask = np.ma.masked_invalid(X_s).mask
        X_s[np.isnan(X_s)] = 0
        
        mask_shared = to_shared(mask, c_bool)
        X_shared = to_shared(X)
        X_s_shared = to_shared(X_s)
        
        pool = mp.Pool(initializer=knn_init, 
                       initargs=(self._k, mask_shared, X_shared, X_s_shared))
        pool.map(knn_work, range(X.shape[0]))
        return from_shared(X_shared)

