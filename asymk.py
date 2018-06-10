#!/bin/python
# This code is an implementation of an approximation algorithm for solving
# Asymmetric k-center problem.
# Given a set V of n points and the distances between each pair
# the k-center problem asks to choose a subset C (subset V) of size k
# that minimizes the maximum over all points of the distance from C
# to the point.
# Asymmetric k-center is a generalization of the k-center where the distances
# are asymmetric.

import numpy as np
import math
from numpy import linalg
from scipy.optimize import linprog

def onehot(length, i):
    ret = np.zeros(length, dtype=int)
    ret[i] = 1
    return ret

def DeriveGraph(D, R):
    """
    Return Gr = (V, Er) where Er = {(u, r) : d(u, v) <= R}
    """
    Gbool = D <= R
    G = Gbool.astype(int)
    G = G - np.identity(D.shape[0])
    print("G = ", G)
    return G

def Gamma(Gr, P, R, isPlus=True):
    """
    if isPlus=True,  Return a set of nodes reachable FROM a node in P within R steps.
    if isPlus=False, Return a set of nodes reachable TO a node in P within R steps.
    P: Set of nodes.
    """
    M = Gr.shape[0]
    GrI = Gr + np.identity(M)
    GrR = np.linalg.matrix_power(GrI, R)

    if isPlus:
        nPaths = np.matmul(P, GrR)
    else:
        nPaths = np.matmul(GrR, P.transpose())
    reachable = nPaths >= 1
    ret = reachable.astype(int)

    # print('M=',M)
    # print('GrI=',GrI)
    # print('GrR=',GrR)
    # print('nPaths=',nPaths)
    # print('reachable=',reachable)

    return ret.flatten()

def CCV(Gr, A):
    """
    Return a center capturning vertex (CCV) in A.
    A node v is CCV iff GammaMinus(v) is a subset of GammaPlus(v).
    """
    for i, val in enumerate(A):
        if val == 1:
            v = onehot(len(A), i)
            GammaPlus = Gamma(Gr, v, 1, True)
            GammaMinus = Gamma(Gr, v, 1, False)
            if np.all(GammaPlus >= GammaMinus):
                # v is a CCV node!
                return v

    # If no CCV node is found, return None.
    return None

def Reduce(Gr):
    """
    """
    M = Gr.shape[0]
    A = np.ones((M,), dtype=int)
    C = np.zeros((M,), dtype=int)

    v = CCV(Gr, A)
    print('v=', v)
    while v is not None:
        # print("Reduce: A=", A)
        # print("Reduce: C=", C)
        C = C + v
        gammaV = Gamma(Gr, onehot(len(C), v), 2)
        for i, val in enumerate(gammaV):
            if val == 1:
                A[i] = 0
        v = CCV(Gr, A)
    # print("Reduce: A=", A)
    # print("Reduce: C=", C)
    gammaC = Gamma(Gr, C, 4)
    for i in range(len(A)):
        if A[i] == 1 and gammaC[i] == 1:
            A[i] = 0
            
    # print("Reduce: A=", A)
    # print("Reduce: C=", C)
    return C, A

def GetGhat(G3r, Gr, C):
    """
    Add G3r an edges of {(u, v): u in C, v in GammaPlus(Gr, C, 4)}
    """
    vs = Gamma(Gr, C, 4, True)
    ghat = np.copy(G3r)
    for i in range(len(C)):
        for j in range(len(vs)):
            if (C[i] == 1) and (vs[j] == 1) and (i is not j):
                ghat[i][j] = 1
    return ghat

def LP(Ghat, A):
    """
    Linear programming to solve Fractional set cover problem:
    min  y(V)
    s.t. y(GammaMinus(v)) >= 1 for all v in A
         y >= 0

    min  sum_{v in V} (y[v])
    s.t. for all v in A: sum_{v' in G-[v]} >= 1
         for all v in V: y[v] >= 0
    Return the assignment of y for the nodes in V. 
    """
    nV = Ghat.shape[0]

    # min y(V)
    # c: Cost function (array of ones)
    c = np.ones_like(A)
    
    # Constraint 1. 
    # y(GammaMinus(v)) >= 1 for all v in A
    constA = A.copy()
    varA = Ghat.transpose()

    # Constraint 2.
    # y >= 0
    constB = np.zeros(nV, dtype=int)
    varB = np.identity(nV)

    
    consts = np.concatenate([-constA, -constB])
    var = np.concatenate([-varA, -varB])

    res = linprog(c, A_ub=var, b_ub=consts)

    print("Solution of LP = ", res.x)
    print("fun = ", res.fun)
    
    return res.x, res.fun

def Vgeqi(G, C, i):
    ret = np.ones_like(C)
    gammaplus = Gamma(G, C, i-1, True)
    ret = ret - gammaplus
    return ret

def Augment(Ghat, A, C, y, p):
    """
    Augment greedly finds a set of centers to add to the already found centers C,
    using the LP-relaxed solution y.
    The algorithm here is based on ExpandingFront (Fig. 4)
    """
    i = 0
    if np.dot(y, A) < 1:
        return C
    while True:
        jmax = np.int(math.ceil(3.0 / 4.0 * p / math.pow(2, i)))
        print('jmax=', jmax)
        for j in range(1, jmax):
            if np.dot(y, A) < 1:
                return C
            else:
                Vip1 = Vgeqi(Ghat, C, i+1)
                maxy = 0
                maxv = -1
                for i, v in enumerate(Vip1):
                    if v == 1:
                        vonehot = onehot(len(A), i)
                        gammav = Gamma(Ghat, vonehot, 1, True)
                        gammavANDA = np.minimum(gammav, A)
                        y = np.dot(y, gammavANDA)
                        if y > maxy:
                            maxy = y
                            maxv = i
                v = maxv
                C[v] = 1 # C <- C + v
                A = Vgeqi(Ghat, C, i+2)
        A = Vgeqi(Ghat, C, i+3)
        i += 1
    return C

def AKC(D, k, R):
    """
    Decision problem of asymmetric k-center.
    V: Nodes (V is implicitly represented in D.)
    D: matrix specifying a distance function V x V -> R.
    k: number of centers.
    R: maximum radius.
    """
    print("################")
    print("AKC: DeriveGraph(D, R)")
    Gr = DeriveGraph(D, R)
    print("Gr = ", Gr)
    
    print("################")
    print("AKC: Reduce(Gr)")
    C, A = Reduce(Gr)
    p = 2.0 / 3.0 * (k - np.sum(C))
    print("C = ", C)
    print("A = ", A)
    print("p = ", p)
    print("################")
    print("AKC: DeriveGraph(D, 3 * R)")
    G3r = DeriveGraph(D, 3 * R)
    print("G3r = ", G3r)
    print("################")
    print("AKC: GetGhat")
    Ghat = GetGhat(G3r, Gr, C)
    print("Ghat = ", Ghat)
    
    print("################")
    print("AKC: LP(Ghat, A)")
    y, yV = LP(Ghat, A)
    print("y = ", y)
    print("y(V) = ", yV)
    if yV > p:
        # R < R*
        print(R, " < R*")
        return None
    else:
        p = yV
        print("################")
        print("AKC: Augment(Ghat, A, C, y, p)")
        C = Augment(Ghat, A, C, y, p)
        return C

def AKC_OPT(D, k):
    print("TODO: AKC_OPT is implemented by incrementing the R one by one.")
    print("This is not efficient: binary search over R.")
    R = 0
    centers = AKC(D, k, R)
    while centers is None:
        R += 1
        centers = AKC(D, k, R)
        
    return AKC(D, k, R), R

def GetRadius(D, C):
    nV = D.shape[0]

    maxd = -1
    for i in range(nV):
        mind = 10000
        for c in range(nV):
            if C[c] == 1:
                dic = D[c][i]
                if dic < mind:
                    mind = dic
        if mind > maxd:
            maxd = mind

    return maxd

if __name__ == "__main__":    
    D = np.array([[0, 1, 2], [100, 0, 1], [100, 100, 0]])
    k = 1
    # C, R = AKC_OPT(D, k)
    C = AKC(D, 1, 1)
    R = GetRadius(D, C)
    
    print("Centers=", C)
    print("Radius=", R)
