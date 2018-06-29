import numpy as np
import functions as fn

# 2x2 Pauli matrices
spin1_x = 1/2 * np.array([[0, 1],
                          [1, 0]])

spin1_y = 1/2 * np.array([[0, -1j],
                          [1j, 0]])

spin1_z = 1/2 * np.array([[1, 0],
                          [0, -1]])

# 4x4 Matrices for S operator
spin2_s_x = fn.kron_a_n(spin1_x, 2)
spin2_s_y = fn.kron_a_n(spin1_y, 2)
spin2_s_z = fn.kron_a_n(spin1_z, 2)
spin2_s_p = spin2_s_x + 1j * spin2_s_y
spin2_s_m = spin2_s_x - 1j * spin2_s_y

# 4x4 Matrices for I operator
spin2_i_x = fn.kron_n_a(2, spin1_x)
spin2_i_y = fn.kron_n_a(2, spin1_y)
spin2_i_z = fn.kron_n_a(2, spin1_z)
spin2_i_p = spin2_i_x + 1j * spin2_i_y
spin2_i_m = spin2_i_x - 1j * spin2_i_y
spin2_all = [spin2_s_z, spin2_s_p, spin2_s_m,
             spin2_i_z, spin2_i_p, spin2_i_m]
