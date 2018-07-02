import multiprocessing
import time
import spin_matrices as sp
import numpy as np


def runSimulation(params):
    """This is the main processing function. It will contain whatever
    code should be run on multiple processors.

    """
    param1, param2 = params
    processedData = np.zeros((int(1E4), 64, 64))
    for ctr in range(int(1E4)):
        for test in range(int(1E1)):
            processedData[ctr] = np.matmul(np.kron(param1, param2), np.kron(param2, param1))
            processedData[ctr] = np.matmul(np.kron(param1, param2), np.kron(param2, param1))
            processedData[ctr] = np.matmul(np.kron(param1, param2), np.kron(param2, param1))


    return processedData


def runSimulation2(param1, param2):
    """This is the main processing function. It will contain whatever
    code should be run on multiple processors.

    """
    processedData = np.zeros((int(1E4), 64, 64))
    for ctr in range(int(1E4)):
        for test in range(int(1E1)):
            processedData[ctr] = np.matmul(np.kron(param1, param2), np.kron(param2, param1))
            processedData[ctr] = np.matmul(np.kron(param1, param2), np.kron(param2, param1))
            processedData[ctr] = np.matmul(np.kron(param1, param2), np.kron(param2, param1))

    return processedData

if __name__ == '__main__':
    # Define the parameters to test
    param1 = sp.spin3_s1_x
    param2 = sp.spin3_s2_x

    params = zip(param1, param2)

    pool = multiprocessing.Pool()

    # # Parallel map
    tic = time.time()
    results2 = pool.map(runSimulation, params)
    toc = time.time()
    print('time', toc - tic)

    # Serial map
    # tic2 = time.time()
    # results2 = map(runSimulation, params)
    # toc2 = time.time()
    # print('time', toc2 - tic2)

    # Serial map
    tic3 = time.time()
    results3 = runSimulation2(param1, param2)
    toc3 = time.time()

    print('time', toc3-tic3)
    print('results2', results2)
    print('results3', results3)
    print('results3-results3', results3-results3)



