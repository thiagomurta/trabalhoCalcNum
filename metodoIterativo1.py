import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.io
from scipy.linalg import lu
import scipy.sparse.linalg

def jacobi(A, b, tol, nmaxiter):
    nColunas = A.shape[0]
    iter = 0
    i = 0
    j = 0
    n = len(b)
    x1 = []
    while(iter < nmaxiter):
        iter += 1
        for i in list(range(0,n+1,1)):
            soma = 0
            for j in list(range(0, n+1, 1)):
                if((i-1) != (j-1)):
                    soma = soma + A[i-1][j-1] * x1[j-1]
            x1[i-1] = (b[i-1] - soma)/A[i-1][i-1]

        erro = np.linalg.norm(x1[i] - x1[i-1])/np.linalg.norm(x1[i-1])
        if(erro < tol):
            return x1, erro, iter


def sor(A, b, tol, nmaxiter, w):
    n = len(b)
    j = 0
    i = 0
    k = 0
    x2 = []
    while j < nmaxiter:
        x_novo = np.copy(x2)
        for i in range(n):
            temp = sum(A[i, j]) * x_novo[j]
            for j in range(n): 
                if j != i:
                    x2[k] = ((1 - w) * x2[k-1]) + (w/A[i,i]) * (b[i] - temp)

            erro = np.linalg.norm(x2 - x_novo)/np.linalg.norm(x_novo)
            if(erro < tol):
                er2 = erro
                return x2, er2, j
        x2 = x_novo 
        j+=1

def fatora(A, w):
    w = 0

def diagonal_dominante(A):
    i = 0
    j = 0
    n = A.shape[0]
    while(i < n):
        if(abs(A[i,i]) > somaLinhaMatriz(i, n, A)):
            i += 1
        else:
            return False
    return True

def somaLinhaMatriz(i, nColunas, A):
    j = 0
    soma = 0
    while(j < nColunas):
        if j != i:
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