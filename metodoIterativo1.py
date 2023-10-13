import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.io
from scipy.linalg import lu
import scipy.sparse.linalg

def jacobi(A, b, tol, nmaxiter):
    nColunas = A.shape[0]
    j = 0
    i = 0
    k = 0
    x1 = []
    while(j < nmaxiter):
        for k in range(nColunas):
            for i in range(nColunas):
                x1[k] = (b[i] - somaLinhaMatrizJacobi(i, nColunas, A, x1) - somaLinhaMatrizJacobi(i+1, nColunas, A, x1))/A[i,i]
        if(abs(x1[k] - x1[k-1]) < tol):
            er = abs(x1[k] - x1[k-1])
            return x1, er, j
    j += 1


def sor(A, b, tol, nmaxiter, w):
    w = 0
    nColunas = A.shape[0]
    j = 0
    i = 0
    k = 0
    x2 = []
    while(j < nmaxiter):
        for k in range(nColunas):
            for i in range(nColunas):
                x2[k] = ((1 - w) * x2[k-1]) + (w/A[i,i]) * (b[i] - somaLinhaMatrizJacobi(i, nColunas, A, x2) - somaLinhaMatrizJacobi(i+1, nColunas, A, x2))
        if(abs(x2[k] - x2[k-1]) < tol):
            er2 = abs(x2[k] - x2[k-1])
            return x2, er2, j
    j+=1

def fatora(A, w):
    w = 0

def diagonal_dominante(A):
    i = 0
    j = 0
    n = A.shape[0]
    m = n
    soma = 0
    while(i < m):
        while(j < m):
            if(abs(A[i,i]) > somaLinhaMatriz(i, n,A)):
                j += 1
            else:
                return False
            m -= 1
            j+=1
        i+=1
    return True

def somaLinhaMatrizJacobi(i, nColunas, A, x):
    j = 0
    soma = 0
    while(j < nColunas):
        soma += A[i,j] * x[j]
        j+=1
    return soma

def somaLinhaMatriz(i, nColunas, A):
    j = 0
    soma = 0
    while(j < nColunas):
        soma += A[i,j]
        j+=1
    return soma

## 01
S = sc.io.mmread('bcsstk02.mtx')
A = S.todense(order='C')
n = A.shape[0]
S = sc.sparse.csr_matrix(A)
print(A)

## 02
b = A.dot(np.ones((n,1)))

## 03
print("A matriz é diagonal Dominante?", diagonal_dominante(A))

## 04
tol = 0.0001
nmaxiter = 10
w = 1.0
x1 = []
x2 = []

x1, er1, iter = jacobi(A, b, tol, nmaxiter)
x2, er2, iter = sor(A, b, tol, nmaxiter, w)
MJ, MS, MSOR = fatora(A, w)

## 05
t, V = np.linalg.eig(A)
max(abs(t))

## 06
y1 = np.log(er1)
plt.plot(x1, y1)

y2 = np.log(er2)
plt.plot(x2, y2)

## 07