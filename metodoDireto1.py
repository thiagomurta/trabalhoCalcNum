import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.io
from scipy.linalg import lu
import scipy.sparse.linalg

## 01
S = sc.io.mmread('bcsstk05.mtx')
A = S.todense(order='C')
n = A.shape[0]
S = scipy.sparse.csr_matrix(A)

## 02 e 03
P, L, U = sc.linalg.lu(A)

## 04
taxaPreenchimento = 100 - ( np.count_nonzero(A) / ( np.count_nonzero(L) + np.count_nonzero(U) ) ) * 100
print("Taxa de Preenchimento: ",taxaPreenchimento)

## 05
u = np.ones((n, 1))
b = A.dot(u)

x = np.linalg.solve(A,b)

## 06
u = np.ones((n, 1))
b = S.dot(u)

y = sc.sparse.linalg.spsolve(S, b)

## 07
distanciaRelativax = np.linalg.norm(u-x)
print("Distância Relativa Linear denso: ", distanciaRelativax)

distanciaRelativay = np.linalg.norm(u-y)
print("Distância Relativa Linear esparso: ", distanciaRelativay)

## 08
alfaA = A - P @ L @ U
#print(alfaA)
t = np.linalg.norm(alfaA, np.inf)
m = np.linalg.norm(A, np.inf)
resultado = t/m
print("Distância Relativa: ",resultado)

## 09
r = np.linalg.norm(u - x)
print("Norma do resíduo: ", r)

## 10
K = np.linalg.cond(A)
print("Número de Condicionamento: ", K)