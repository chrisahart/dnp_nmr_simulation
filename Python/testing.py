import time
import spin_matrices as sp
import numpy as np
from scipy import linalg as la
import f2py_dynamics as fortran

A = sp.spin3_i_x
B = sp.spin3_i_y
C = int(1E6)

# A=([[1, 2, 3, 4],
#     [5, 6, 7, 8],
#     [9, 10, 11, 12],
#     [13, 14, 15, 16]])
#
# B=np.array([[1, 5, 9, 13],
#     [2, 6, 10, 14],
#     [3, 7, 11, 15],
#     [4, 8, 12, 16]])
#
# Ix = np.kron(A, B)
# Iy = np.kron(B, A)

start = time.time()
fortran.f2py_dynamics.calculate_hamiltonian()
end = time.time() - start
print('time taken', end)


# for count in range(C):
#     temp = np.matmul(A, B)
#     # temp = np.kron(A, B)
#     # temp = la.expm(-1j * sp.spin3_i_x * int(1E4))
# end = time.time() - start
# print('time taken', end)

