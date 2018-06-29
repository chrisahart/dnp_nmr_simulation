from multiprocessing import Pool
import spin_matrices as sm
import numpy as np
import parameters as param
import time
import functions as fn

# def process_image(name):
#
#     b = a^2
#
#     return
#
#
# if __name__ == '__main__':
#     pool = Pool()                         # Create a multiprocessing Pool
#
#     a = np.linspace(0, 100, num=101)
#
#     pool.map(process_image, data_inputs(a))  # process data_inputs iterable with pool

# A = sm.spin2_s_x
# N = 2 ** 2
#
# start = time.time()
#
# for count1 in range(0, int(param.time_step_num)):
#     C = np.kron(A, np.eye(N))
#
# print('time taken', time.time() - start)
# print('C \n', C.shape)
#
# start = time.time()
#
# for count2 in range(0, int(param.time_step_num)):
#     D = fn.kron_a_n(A, N)
#
# print('time taken', time.time() - start)
# print('D \n', D.shape)
#
# print('C-D \n', C-D)