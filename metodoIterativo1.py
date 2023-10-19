import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import scipy.io
from scipy.linalg import lu
import scipy.sparse.linalg

# Função Jacobi
def jacobi(A, b, tol, nmaxiter):
    nColunas = A.shape[0]
    iter = 0
    n = len(b)
    errovet = []
    x1 = np.zeros((n, n))
    for k in range(n):
        soma = 0
        x_old = x1.copy()
        for i in range(nColunas):
            iter += 1
            erro = 0
            soma += (somaLinhaMatrizInferior(i, A) * x_old[k] ) + (somaLinhaMatrizSuperior(i, A) * x_old[k] )
            x1[k] = (b[k] - soma)/A[i,i]
            erro = np.linalg.norm(x1 - x_old)
            errovet.append(erro)
            if(np.any(erro < tol) or iter >= nmaxiter):
                return x1, errovet, iter

# Função SOR
def sor(A, b, tol, nmaxiter, w):
    nColunas = A.shape[0]
    n = len(b)
    iter = 0
    errovet = []
    x1 = np.zeros((n, n))
    for k in range(n):
        soma = 0
        x_novo = x1.copy()
        for i in range(nColunas):
            iter += 1
            erro = 0
            soma = ( somaLinhaMatrizInferior(i, A) * x_novo[i] ) + (somaLinhaMatrizSuperior(i, A) * x_novo[i-1])
            x1[k] = ((1 - w) * x1[k-1]) + (w/A[i,i]) * (b[k] - soma)
            erro = np.linalg.norm(x1 - x_novo)
            errovet.append(erro)
            if(np.any(erro < tol) or iter >= nmaxiter):
                return x1, errovet, iter

# Função fatora
def fatora(A, w):
    n = A.shape[0]
    D = np.diag(np.diag(A))
    L = np.tril(A, k=-1)
    U = np.triu(A, k=1)
    # Matrizes
    MJ = np.linalg.inv(D).dot(L+U)
    MS = np.linalg.inv(D - L).dot(U)
    MSOR = np.linalg.inv(D + w * L).dot((1 - w) * D - w * U)
    return MJ, MS, MSOR

# Função pra verificar se a matriz é diagonal dominante
def diagonal_dominante(A):
    i = 0
    n = A.shape[0]
    while(i < n):
        if(abs(A[i,i]) > somaLinhaMatriz(i, A)):
            i += 1
        else:
            return False
    return True

# Funções auxiliares
def somaLinhaMatriz(i, A):
    nColunas = A.shape[0]
    j = 0
    soma = 0
    for j in range(nColunas):
        if j != i:
            soma += A[i,j]
    return soma

def somaLinhaMatrizInferior(i, A):
    nColunas = A.shape[0]
    j = 0
    soma = 0
    for j in range(nColunas):
        if j != i and j < i:
            soma += A[i,j]
    return soma

def somaLinhaMatrizSuperior(i, A):
    nColunas = A.shape[0]
    j = 0
    soma = 0
    for j in range(nColunas):
        if j != i and j > i:
            soma += A[i,j]
    return soma

## 01
S = sc.io.mmread('bcsstk05.mtx')
A = S.todense(order='C')
n = A.shape[0]
S = sc.sparse.csr_matrix(A)
print(A)

## 02
b = A.dot(np.ones((n,1)))

## 03
print("A matriz é diagonal Dominante?", diagonal_dominante(A))


## 04
w = 1.25
MJ, MS, MSOR = fatora(A, w)
raioJacobi, T = np.linalg.eig(MJ)
raioSeidel, S = np.linalg.eig(MS)
raioSOR, M = np.linalg.eig(MSOR)
print("Raio espectral Jacobi: ", max(abs(raioJacobi)))
print("Raio espectral Seidel: ", max(abs(raioSeidel)))
print("Raio espectral SOR: ", max(abs(raioSOR)))

## 05
tol = 1e-6
nmaxiter = 100
x1 = []
x2 = []
er1 = []
er2 = []
if np.any(raioSOR < 1.0):
    x1, er1, iter1 = jacobi(A, b, tol, nmaxiter)
    x2, er2, iter2 = sor(A, b, tol, nmaxiter, w)
    print("Interações Jacobi:", iter1)
    print("Interações SOR:", iter2)

## 06
# Eixo x e y do Jacobiano
xJac = np.linspace(1, iter1, num=iter1,endpoint=True)
y1 = np.log(er1)

# Eixo x e y do Gauss Seidel
xGS = np.linspace(1, iter2, num=iter2,endpoint=True)
y2 = np.log(er2)

# Criando o gráfico
plt.plot(xJac, y1)
plt.plot(xGS, y2)
plt.grid(True)
plt.xlabel('Iterações(i)')
plt.ylabel('Log(erro)')
plt.show()

## 07