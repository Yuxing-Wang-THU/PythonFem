import numpy as np

A = np.array([[5, 2, 1],
     [-1, 4, 2],
     [2, -3, 10]])
bb = np.array([[-12], [20], [3]])

def JacoBinSolver(A, b, iter=100, error=0.001):
    # Formulate DLU
    m, n = np.shape(A)
    D = np.mat(np.zeros((m, n)))
    L = np.mat(np.zeros((m, n)))
    U = np.mat(np.zeros((m, n)))
    for i in range(m):
        for j in range(n):
            if i == j:
                D[i, j] = A[i, j]
            if i < j:
                L[i, j] = -A[i, j]
            if i > j:
                U[i, j] = -A[i, j]
    
    b = np.reshape(b,(-1,1))
    # Initial value
    x0 = np.array(np.zeros((m, 1)))
    xk = np.array(np.zeros((m, 1)))

    B = np.dot(D.I, (L + U))
    f = np.dot(D.I, b)

    iter_time = 1
    xk = np.dot(B, x0) + f
    while(np.linalg.norm((xk - x0)) >= error):
        iter_time += 1
        x0 = xk
        xk = np.dot(B, xk) + f
        if iter_time > iter:
            break
    result = np.squeeze(xk.A)
    return result 


if __name__ == "__main__":
    aa = JacoBinSolver(A,bb)
    print(aa)