import time
import spin_matrices as sp
import numpy as np
from scipy import linalg as la

A = sp.spin3_i_x
B = sp.spin3_i_y

start = time.time()

for count in range(int(1E5)):
    # temp = np.matmul(A, B)
    temp = np.kron(A, B)
    # temp = la.expm(-1j * sp.spin3_i_x * int(1E4))

print('time taken', (time.time() - start))
