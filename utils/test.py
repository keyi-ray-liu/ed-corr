import numpy as np
from numpy.linalg import eig


N = 100
# M = np.random.random((N, N))

# M = (M + np.conj(np.transpose(M)))
# w, v = eig(M)

# print(np.allclose( M, v @ np.diag(w) @ np.conj(np.transpose(v))))

a = np.random.random(N)
A = np.diag(a)
B = np.random.random((N, N))

print( A @ B - B @ A)