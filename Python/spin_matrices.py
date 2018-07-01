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

# 8x8 Matrix operators for S1 operator
spin3_s1_x = np.kron(np.kron(spin1_x, np.identity(2)), np.identity(2))
spin3_s1_y = np.kron(np.kron(spin1_y, np.identity(2)), np.identity(2))
spin3_s1_z = np.kron(np.kron(spin1_z, np.identity(2)), np.identity(2))
spin3_s1_p = spin3_s1_x + 1j * spin3_s1_y
spin3_s1_m = spin3_s1_x - 1j * spin3_s1_y


# 8x8 Matrix operators for S2 operator
spin3_s2_x = np.kron(np.identity(2), np.kron(spin1_x, np.identity(2)))
spin3_s2_y = np.kron(np.identity(2), np.kron(spin1_y, np.identity(2)))
spin3_s2_z = np.kron(np.identity(2), np.kron(spin1_z, np.identity(2)))
spin3_s2_p = spin3_s2_x + 1j * spin3_s2_y
spin3_s2_m = spin3_s2_x - 1j * spin3_s2_y

# 8x8 Matrix operators for I operator
spin3_i_x = np.kron(np.identity(2), np.kron(np.identity(2), spin1_x))
spin3_i_y = np.kron(np.identity(2), np.kron(np.identity(2), spin1_y))
spin3_i_z = np.kron(np.identity(2), np.kron(np.identity(2), spin1_z))
spin3_i_p = spin3_i_x + 1j * spin3_i_y
spin3_i_m = spin3_i_x - 1j * spin3_i_y
spin3_all = [spin3_s1_z, spin3_s1_p, spin3_s1_m,
             spin3_s2_z, spin3_s2_p, spin3_s2_m,
             spin3_i_z, spin3_i_p, spin3_i_m]
