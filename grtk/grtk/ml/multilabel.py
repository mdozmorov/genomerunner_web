# -*- coding: utf-8 -*-
import numpy
import math
from math import sqrt
from numpy import dot, fabs, maximum, log, sign, array, concatenate, zeros
from numpy.linalg import norm
import sklearn.datasets
import psycopg2
from pandas import SparseDataFrame

def logistic(x):
    return numpy.power(1 + numpy.exp(-1 * x), -1)

def l21norm(X):
    r = 0
    for i in range(X.shape[0]):
        r += norm(X[i,:], 2)
    return r

class MultilabelLR(object):
    """
    This is an implementation of an approach to multilabel logistic regression
    with dense conditional dependency networks, as described in the paper:

    Yuhong Guo and Wei Xue. Probabilistic multi-label classification with
    sparse feature learning. IJCAI 2013.
    
    http://www.cis.temple.edu/~yuhong/research/papers/ijcai13b.pdf
    """
    def __init__(self, λ1=1e-1, λ2=1e-1, γ=1e-4, μ=1e-5, η=2, L0=1):
        """
        - L0 : Initial guess for Lipschitz constant
        - η  : search rate for the Lipschitz constant (must be > 1)
        - λ1, λ2 : regularization parameters for the coefficients on
                   X and Y, respectively
        """
        assert(η > 1)
        #assert(λ1 > 0)
        #assert(λ2 > 0)
        assert(γ > 0)
        assert(μ > 0)
        assert(L0 > 0)

        self.λ1 = λ1
        self.λ2 = λ2
        self.γ = γ
        self.μ = μ
        self.η = η
        self.L0 = L0
    
    def decompose(self, θ):
        W = θ[:self.D,:]
        V = θ[self.D:-1,:]
        b = θ[-1,:]
        return W,V,b

    def fit(self, X, Y):
        """
        Input Variables:
        - X : a 2D array of instances (rows) and features (columns)
        - Y : a 2D array of instances (rows) and labels (columns)

        Internal Variables:
        - K    : # of labels
        - D    : # of features
        - W    : [D x K] shaped matrix of coefficients for training features 
                 (within X)
        - V    : [K-1 x K] shaped matrix of coefficients reflecting inter-label
                 correlations
        - b    : K length vector of bias terms
        - θ    : the full coefficient matrix resulting from the vertical 
                 concatenation of W, V, and b
        """
        L = self.L0
        N = X.shape[0]
        α = 1
        self.D = D = X.shape[1]
        self.K = K = Y.shape[1]
        z = numpy.zeros((D+K,K))

        def f(θ):
            """
            The smooth, convex part of the objective function F,
            consisting of the log-likelihood of the logistic function and
            regularization terms.
            """
            W,V,b = self.decompose(θ)
            f = 0
            for k in range(K):
                w = W[:,k]
                v = V[:,k]
                y = Y[:,k]
                z = y * dot(X, w) + dot(Y[:,:k], v[:k]) \
                    + dot(Y[:,(k+1):], v[k:]) + b[k]
                f += log(logistic(z)).sum()
            f *= -1
            f += (self.λ1 / 2) * pow(norm(W, 'fro'), 2)
            f += (self.λ2 / 2) * pow(norm(V, 'fro'), 2)
            return f

        def Δf(θ):
            "The gradient of f at point theta."
            W, V, b = self.decompose(θ)
            gradient = numpy.zeros(θ.shape)
            for k in range(K):
                z = dot(X, W[:,k]) + dot(Y[:,:k], V[:k,k]) + \
                    dot(Y[:,(k+1):], V[k:,k]) + b[k]
                l = logistic(-1 * Y[:,k] * z)
                for t in range(D):
                    gradient[t,k] = (Y[:,k] * X[:,t] * l).sum()
                for t in range(D, D+k):
                    gradient[t,k] = (Y[:,k] * Y[:,t-D] * l).sum()
                for t in range(D+k+1, D+K):
                    gradient[t-1,k] = (Y[:,k] * Y[:,t-D] * l).sum()
                gradient[-1,k] = (Y[:,k] * l).sum()

            gradient[:] *= -1
            gradient[:self.D,:] -= self.λ1 * W
            gradient[self.D:-1,:] -= self.λ2 * V
            return gradient 
 
        def g(θ):
            """
            The non-smooth part of the objective function, consisting
            of sparsity-inducing regularization terms.
            """
            W,V,b = self.decompose(θ)
            g = 0
            for k in range(K):
                g += norm(V[:,k], 1)
            return (g * self.γ) + self.μ * l21norm(W)

        def F(θ):
            "The objective function."
            return f(θ) + g(θ)

        def qL(θ, θt):
            #print(dot((θ - θt).flat, Δf(θt).flat))
            return (f(θt) + dot((θ - θt).flat, Δf(θt).flat) + \
                (L / 2) * pow(norm((θ - θt).flat, 2), 2) + g(θ))

        def pL(θ):
            "Compute the next update given the current Lipschitz constant."
            theta_hat = θ - Δf(θ) / L
            W_hat, V_hat, b_hat = self.decompose(theta_hat)
            for j in range(D):
                W_hat[j,:] *= maximum(1 - (self.μ / (L * norm(W_hat[j,:], 2))), 0)
            for k in range(K):
                V_hat[:,k] = sign(V_hat[:,k]) \
                             * maximum(0, fabs(V_hat[:,k]) - self.γ / L)
            return theta_hat

        """
        pLz = pL(z)
        while True: #(fabs(pLz - z).sum() > 1e-6):
            pLz = pL(z)
            print("L =", L, "F(z) =", F(z), "F(pL(z)) =", F(pLz))
            if F(pLz) >= F(z):
                L = self.η * L
            else:
                z = pLz
            if L > 1000000:
                break
            
        self.θ = pLz
        print("θ sparsity % =", (pLz==0).sum() / (len(pLz.flat)))
                
        """
        for i in range(1000):
            print("\n*** Iteration #", i+1)
            while True:
                pLz = pL(z)
                if F(pLz) <= qL(pLz, z):
                    break
                print("L =", L, "F(z) =", F(z), 
                      "F(pL(z)) =", F(pLz), "qL =", qL(pLz, z))
                L = self.η * L
            α_next = (1 + sqrt(1 + 4 * pow(α, 2))) / 2
            z_next = pLz + ((α - 1) / α_next) * (pLz - z)
            dx_avg = fabs(z_next - z).mean()

            print("F =", f(pLz))
            print("g =", g(pLz))
            print("α =", α)
            print("L =", L)
            print("Mean delta =", dx_avg)

            α = α_next
            z = z_next
            if L > 1e7:
                break
        
        self.θ = pLz
        print("θ sparsity % =", (pLz==0).sum() / (len(pLz.flat)))

    def _p(self, X, Y, θ, P=None, order=None):
        W,V,b = self.decompose(θ)
        if P is None:
            P = zeros(Y.shape)
        if order is None:
            order = range(self.K)
        for k in order:
            w = W[:,k]
            v = V[:,k]
            y = Y[:,k]
            z = y * dot(X, w) + dot(Y[:,:k], v[:k]) \
                + dot(Y[:,(k+1):], v[k:]) + b[k]
            P[:,k] = logistic(z)
        return P

    def predict_proba(self, X):
        burn_in = 100
        samples = 100
        skip = 5

        W, V, b = self.decompose(self.θ)
        N = X.shape[0]
        Y = numpy.random.choice([-1,1], size=N * self.K).reshape((N,self.K))
        P = zeros((Y.shape))

        order = array(range(self.K))
        numpy.random.shuffle(order)

        Y_count = zeros((N,self.K))
        for i in range(burn_in + samples * skip):
            P = self._p(X, Y, self.θ, P=P, order=order)
            Y[:] = 1
            Y.flat[numpy.random.random(N * self.K) > P.flat] = -1
            if (i > burn_in) and (i % skip == 0):
                Y_count[Y==1] += 1
        return Y_count / samples
    
    def predict(self, X):
        P = self.predict_proba(X)
        for i in range(P.shape[0]):
            for j in range(P.shape[1]):
                P[i,j] = 1 if P[i,j] > 0.5 else -1
        return P


from sklearn.datasets import make_multilabel_classification
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelBinarizer
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
    
def test():
    X,Y_list = make_multilabel_classification()
    Y = LabelBinarizer().fit_transform(Y_list)
    Y[Y==0] = -1
    clf = OneVsRestClassifier(LinearSVC())
    #clf = MultilabelLR(L0=1, λ1=0.1, λ2=0.1, γ=0.1, μ=0.1)
    clf.fit(X,Y)
    Y_hat = clf.predict(X)

    print(roc_auc_score(Y.flat, Y_hat.flat))

test()

#ix = ((Ts == 1).sum(axis=0) > 500)
#clf = MultilabelLR(L0=1, λ1=0.1, λ2=0.1, γ=0.1, μ=0.1)

#clf = MultilabelLR()
#clf.fit(Xs[:1000,:100],Ts[:1000,ix])
#T_ss = Ts[:100,ix]
#T_hat = clf.predict_proba(Xs[:100,:100])
#import sklearn.metrics
#auc = sklearn.metrics.roc_auc_score(T_ss[:,5], T_hat[:,5])
#print(auc)

#clf = LogisticRegression()
#clf.fit(Xs[:,:100], Ts[:,ix][:,0])
#y_hat = clf.predict_proba(Xs[:1000,:100])[:,1]
#print(sklearn.metrics.roc_auc_score(Ts[:1000,ix][:,0], y_hat))
