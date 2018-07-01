import numpy as np
import numba
import time as time


@numba.jit()
def add2_par(x, y):
    c = np.zeros(int(1e6))
    for a in range(0, int(1e6)):
        c[a] = x[a] + 2 * y[a]
    return c


n = np.array(int(1e6))
X = np.ones(n, dtype=np.float)
Y = np.ones(n, dtype=np.float)

t0 = time.time()
#add2_par(X, Y, out=X)
c = add2_par(X, Y)
t1 = time.time()
print('c', c)
print('time', t1 - t0)
