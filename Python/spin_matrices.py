import numpy as np

# 2x2 Pauli matrices
spin1_x = 1/2 * np.array([[0, 1],
                          [1, 0]])

spin1_y = 1/2 * np.array([[0, -1j],
                          [1j, 0]])

spin1_z = 1/2 * np.array([[1, 0],
                          [0, -1]])

# 4x4 Matrices for S operator
spin2_s_x = np.kron(spin1_x, np.identity(2))
spin2_s_y = np.kron(spin1_y, np.identity(2))
spin2_s_z = np.kron(spin1_z, np.identity(2))
spin2_s_p = spin2_s_x + 1j * spin2_s_y
spin2_s_m = spin2_s_x - 1j * spin2_s_y

# 4x4 Matrices for I operator
spin2_i_x = np.kron(np.identity(2), spin1_x)
spin2_i_y = np.kron(np.identity(2), spin1_y)
spin2_i_z = np.kron(np.identity(2), spin1_z)
spin2_i_p = spin2_i_x + 1j * spin2_i_y
spin2_i_m = spin2_i_x - 1j * spin2_i_y
