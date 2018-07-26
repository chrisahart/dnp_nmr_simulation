import time
import spin_matrices as sp
import numpy as np
from scipy import linalg as la
import f2py_dynamics as fortran
import parameters as param
import functions as fn
import matplotlib.pyplot as plt

# A = sp.spin3_i_x
# B = sp.spin3_i_z
# C = int(1E6)

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
#     temp = np.dot(A, B)
# end = time.time() - start
# print('time taken', end)
# print('temp', temp)

# matrix = sp.spin2_s_x #sp.spin1_x + sp.spin1_y
# origonal = la.expm(-1j * matrix * int(1E4))
# print('la.expm \n', origonal)

# eigvals, eigvectors = np.linalg.eig(sp.spin2_i_z)
# print(eigvectors)

# hamiltonian = np.zeros((100, 2 ** 2, 2 ** 2,), dtype=np.complex)
#
# for count in range(0, 99):
#
#     hamiltonian[count, :] = (count+1) * 2E6 * sp.spin2_s_z + \
#                             1E6 * sp.spin2_i_z #+ \
#                             #2 * 1E4 * np.matmul(sp.spin2_i_x, sp.spin2_s_z)

# spin_x = 1/2 * np.array([[0, 1],  [1, 0]])
# spin_z = 1/2 * np.array([[1, 0], [0, -1]])
#
# #hamiltonian = 2E6 * spin_z + 1E6 * spin_x + 1E6 * np.matmul(spin_x, spin_z)
# hamiltonian = 2E6 * sp.spin2_s_z + 1E6 * sp.spin2_i_z + 1E6 * np.matmul(sp.spin2_i_x, sp.spin2_s_z)
#
# print('hamiltonian', hamiltonian)
#
# eigvals, eigvectors = np.linalg.eig(hamiltonian)
# eigvectors_inv = np.linalg.inv(eigvectors)
#
# print('hamiltonian', hamiltonian)
# print('eigvals', eigvals)
# print('eigvectors', eigvectors)
# print('eigvectors_inv', eigvectors_inv)

# start = time.time()
# fortran.f2py_dynamics.testing()
# end = time.time() - start
# print('time taken', end)

size = 80
mat = np.ones((size, size))

for a in range(0, size):
    for b in range(0, size):
        mat[a, b] = ((a+1) * (b+1))

mat = mat + 1E-9

start = time.time()
for loop in range(0, int(1E4)):
    test = la.expm(-1j * mat * 1E-8)


end = time.time() - start
print('time taken', end)

matlab = [0.91, 0.97, 2.36, 5.45, 8.69]
python = [2.59, 2.89, 9.70, 35.4, 72.7]
fortran = [0.037, 0.12, 4.14, 30.6, 74.9]
openmp = [0.0039, 0.016, 0.52, 3.87, 9.53]
sizes = [4, 8, 30, 60, 80]

fig_bench = plt.figure()
ax_bench = fig_bench.add_subplot(111)
ax_bench.plot(sizes, matlab, 'rx')
ax_bench.plot(sizes, matlab, 'r', label='Matlab')
ax_bench.plot(sizes, python, 'gx')
ax_bench.plot(sizes, python, 'g', label='Python')
ax_bench.plot(sizes, fortran, 'bx')
ax_bench.plot(sizes, fortran, 'b', label='Fortran (1)')
ax_bench.plot(sizes, openmp, 'kx')
ax_bench.plot(sizes, openmp, 'k', label='Fortran (8)')
ax_bench.set_xlim([4, 80])
ax_bench.set_xlabel('Matrix size')
ax_bench.set_ylabel('Time (s)')
ax_bench.legend()
plt.show()

