import time
import spin_matrices as sp
import numpy as np
from scipy import linalg as la
import f2py_dynamics as fortran

A = sp.spin1_x #sp.spin3_i_x
B = sp.spin1_z #sp.spin3_i_y
C = int(1E7)

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
#
# start = time.time()
# for count in range(C):
#     # temp = np.matmul(A, B)
#     #temp = np.kron(A, B)
#     # temp = la.expm(-1j * sp.spin3_i_x * int(1E4))
#     temp = np.matmul(A, B)
# end = time.time() - start
# print('time taken', end)

# matrix = sp.spin2_s_x #sp.spin1_x + sp.spin1_y
# origonal = la.expm(-1j * matrix * int(1E4))
# print('la.expm \n', origonal)

start = time.time()
fortran.f2py_dynamics.testing()
end = time.time() - start
print('time taken', end)
