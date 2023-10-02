import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.io
from scipy.linalg import lu
import scipy.sparse.linalg

def jacobi(A, b, tol, nmaxiter):
    tol = 0

def sor(A, b, tol, nmaxiter, w):
    w = 0

def fatora(A, w):
    w = 0

def diagonal_dominante(A):
    i = 0
    j = 0
    k = 0
    h = 0
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
S = sc.io.mmread('bcsstk02.mtx')
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