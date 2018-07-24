import numpy as np
import functions as fn
import time

# 2x2 Pauli matrices
spin1_x = 1/2 * np.array([[0, 1], [1, 0]])
spin1_y = 1/2 * np.array([[0, -1j], [1j, 0]])
spin1_z = 1/2 * np.array([[1, 0], [0, -1]])

# 4x4 Matrices for S operator
spin2_s_x = fn.kron_rmat_eye(spin1_x, 2)
spin2_s_y = fn.kron_rmat_eye(spin1_y, 2)
spin2_s_z = fn.kron_rmat_eye(spin1_z, 2)
spin2_s_p = np.real(spin2_s_x + 1j * spin2_s_y)
spin2_s_m = np.real(spin2_s_x - 1j * spin2_s_y)

# 4x4 Matrices for I operator
spin2_i_x = fn.kron_eye_rmat(2, spin1_x)
spin2_i_y = fn.kron_eye_rmat(2, spin1_y)
spin2_i_z = fn.kron_eye_rmat(2, spin1_z)
spin2_i_p = np.real(spin2_i_x + 1j * spin2_i_y)
spin2_i_m = np.real(spin2_i_x - 1j * spin2_i_y)
spin2_all = [spin2_s_z, spin2_s_x, spin2_s_p, spin2_s_m,
             spin2_i_z, spin2_i_x, spin2_i_p, spin2_i_m]

# 8x8 Matrix operators for S1 operator
spin3_s1_x = fn.kron_rmat_eye(fn.kron_rmat_eye(spin1_x, 2), 2)
spin3_s1_y = fn.kron_rmat_eye(fn.kron_rmat_eye(spin1_y, 2), 2)
spin3_s1_z = fn.kron_rmat_eye(fn.kron_rmat_eye(spin1_z, 2), 2)
spin3_s1_p = np.real(spin3_s1_x + 1j * spin3_s1_y)
spin3_s1_m = np.real(spin3_s1_x - 1j * spin3_s1_y)

# 8x8 Matrix operators for S2 operator
spin3_s2_x = fn.kron_eye_rmat(2, fn.kron_rmat_eye(spin1_x, 2))
spin3_s2_y = fn.kron_eye_rmat(2, fn.kron_rmat_eye(spin1_y, 2))
spin3_s2_z = fn.kron_eye_rmat(2, fn.kron_rmat_eye(spin1_z, 2))
spin3_s2_p = np.real(spin3_s2_x + 1j * spin3_s2_y)
spin3_s2_m = np.real(spin3_s2_x - 1j * spin3_s2_y)

# 8x8 Matrix operators for I operator
spin3_i_x = fn.kron_eye_rmat(2, fn.kron_eye_rmat(2, spin1_x))
spin3_i_y = fn.kron_eye_rmat(2, fn.kron_eye_rmat(2, spin1_y))
spin3_i_z = fn.kron_eye_rmat(2, fn.kron_eye_rmat(2, spin1_z))
spin3_i_p = np.real(spin3_i_x + 1j * spin3_i_y)
spin3_i_m = np.real(spin3_i_x - 1j * spin3_i_y)
spin3_all = [spin3_s1_z, spin3_s1_x, spin3_s1_p, spin3_s1_m,
             spin3_s2_z, spin3_s2_x, spin3_s2_p, spin3_s2_m,
             spin3_i_z, spin3_i_x, spin3_i_p, spin3_i_m]
