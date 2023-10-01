import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

## 01
S = sc.io.mmread('bcsstk02.mtx')
A = S.todense(order='C')
n = A.shape[0]

## 02 e 03
P, L, U = sc.linalg.lu(A)

## 04
taxaPreenchimento = 100 - ( nnz(A) / ( nnz(L) + nnz(U)) ) * 100

## 05
u = np.ones((n, 1))
b = A.dot(u)

x = np.linalg.solve(A,b)
## 06
u = np.ones((n, 1))
b = S.dot(u)

x = sc.sparse.linalg.spsolve(S,b)

## 07
distanciaRelativa = abs(u - x)

## 08 ??????????

## 09 ??????????
#r = sqrt(b -)

## 10
K = sc.linalg.cond(A)