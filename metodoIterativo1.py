import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.io
from scipy.linalg import lu
import scipy.sparse.linalg

def jacobi(A, b, tol, nmaxiter):
    nColunas = np.shape[A]
    j = 0
    x = []
    while(j < nmaxiter):
        for k in A:
            for i in A:
                x[k] = (b[i] - somaLinhaMatriz(i, nColunas, A) - somaLinhaMatriz(i+1, nColunas, A))/A[i,i]
        if(abs(x[i] - x[i-1]) < tol):
            return x[i]


def sor(A, b, tol, nmaxiter, w):
    w = 0

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

def somaLinhaMatriz(i, nColunas, A):
    j = 0
    soma = 0
    while(j < nColunas):
        if(j != i):
            soma += A[i,j]
        j+=1
    return soma

## 01
S = sc.io.mmread('bcsstk25.mtx')
A = S.todense(order='C')
n = A.shape[0]
print(A)

## 02
b = A.dot(np.ones((n,1)))

## 03
print("A matriz é diagonal Dominante?", diagonal_dominante(A))

## 04


## 05


## 06


## 07