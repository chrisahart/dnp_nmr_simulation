import numpy as np
import parameters as param
import spin_matrices as sp
from scipy import constants as sc
from scipy import linalg as la
import f2py_dynamics as fortran_se
import matplotlib.pyplot as plt
import f2py_cross_effect as fortran_ce
import time


def liouville_propagator(num_spins, energies, eigvectors, eigvectors_inv,
                         microwave_hamiltonian_init, calculate_relaxation_mat, spin_all):
    """ Calculate Liouville space propagator.
     """

    # Pre-allocate propagator
    propagator = np.zeros((int(param.time_step_num), 4 ** num_spins, 4 ** num_spins), dtype=np.complex)

    # Calculate variables for original Liouville space relaxation
    p_e = np.tanh(0.5 * param.electron_frequency * (sc.Planck / (sc.Boltzmann * param.temperature)))
    p_n = np.tanh(0.5 * param.freq_nuclear_1 * (sc.Planck / (sc.Boltzmann * param.temperature)))
    gnp = 0.5 * (1 - p_n) * (1 / (1 * param.t1_nuc))
    gnm = 0.5 * (1 + p_n) * (1 / (1 * param.t1_nuc))
    gep = 0.5 * (1 - p_e) * (1 / (1 * param.t1_elec))
    gem = 0.5 * (1 + p_e) * (1 / (1 * param.t1_elec))

    # Calculate variables for Mance Liouville space relaxation
    boltzmann_elec = param.electron_frequency * (sc.Planck / (sc.Boltzmann * param.temperature))
    boltzmann_nuc = param.freq_nuclear_1 * (sc.Planck / (sc.Boltzmann * param.temperature))

    # spin3_s1_z = spin_all[1]

    for count in range(0, int(param.time_step_num)):

        # Transform microwave Hamiltonian into time dependent basis
        microwave_hamiltonian = np.matmul(eigvectors_inv[count], np.matmul(microwave_hamiltonian_init,
                                                                           eigvectors[count]))

        # Calculate total Hamiltonian
        total_hamiltonian = np.diag(energies[count]) + microwave_hamiltonian

        # start = time.time()
        # for bench in range(0, int(1E6)):
        #     #test = np.matmul(eigvectors_inv[count], np.matmul(spin3_s1_z, eigvectors[count]))
        #     # eigvectors_liouville = np.kron(eigvectors[count], eigvectors[count])
        #     hamiltonian_liouville = kron_rmat_eye(total_hamiltonian, 2 ** num_spins)
        # print('kron_rmat_eye', (time.time() - start))
        #
        # start = time.time()
        # for bench in range(0, int(1E6)):
        #     # test = np.matmul(eigvectors_inv[count], np.matmul(spin3_s1_z, eigvectors[count]))
        #     # eigvectors_liouville = np.kron(eigvectors[count], eigvectors[count])
        #     hamiltonian_liouville = np.kron(total_hamiltonian, np.eye(4))
        # print('np.kron', (time.time() - start))

        # Transform Hilbert space Hamiltonian into Liouville space
        hamiltonian_liouville = kron_rmat_eye(total_hamiltonian, 2 ** num_spins) - \
                                kron_eye_rmat(2 ** num_spins, np.transpose(total_hamiltonian))

        # Calculate time dependent Liouville space relaxation matrix
        relax_mat = calculate_relaxation_mat(eigvectors[count], eigvectors_inv[count], gnp, gnm, gep, gem, spin_all)

        # Calculate time dependent solid effect origonal Liouville space relaxation matrix using F2PY
        # relax_mat = fortran.solid_effect_dynamics.calculate_relaxation_mat(
        #     eigvectors[count], eigvectors_inv[count], np.eye(4), np.eye(16), 16, 4, sp.spin2_s_z, sp.spin2_s_p,
        #     sp.spin2_s_m, sp.spin2_i_z, sp.spin2_i_p, sp.spin2_i_m, param.t2_elec, param.t2_nuc, gnp, gnm, gep, gem)

        # Calculate time dependent cross effect origonal Liouville space relaxation matrix using F2PY
        # relax_mat = fortran_ce.cross_effect_dynamics.calculate_relaxation_mat(
        #     eigvectors[count], eigvectors_inv[count], np.eye(8), np.eye(64), 64, 8, sp.spin3_s1_z, sp.spin3_s1_p,
        #     sp.spin3_s1_m, sp.spin3_s2_z, sp.spin3_s2_p, sp.spin3_s2_m, sp.spin3_i_z, sp.spin3_i_p, sp.spin3_i_m,
        #     param.t2_elec, param.t2_nuc, gnp, gnm, gep, gem)

        # Calculate time dependent cross effect Mance Liouville space relaxation matrix using F2PY
        # relax_mat = fortran_ce.cross_effect_dynamics.calculate_relaxation_mat_mance(
        #     eigvectors[count], eigvectors_inv[count], 64, 8, sp.spin3_s1_x, sp.spin3_s1_z,
        #     sp.spin3_s2_x, sp.spin3_s2_z, sp.spin3_i_x, sp.spin3_i_z, param.t1_elec, param.t1_nuc, param.t2_elec,
        #     param.t2_nuc, boltzmann_elec, boltzmann_nuc)

        # Calculate Liouville space eigenvectors
        eigvectors_liouville = np.kron(eigvectors[count], eigvectors[count])
        eigvectors_inv_liouville = np.kron(eigvectors_inv[count], eigvectors_inv[count])

        # Calculate Liouville space propagator
        liouvillian = hamiltonian_liouville + 1j * relax_mat
        propagator[count, :] = np.matmul(eigvectors_inv_liouville,
                                         np.matmul(la.expm(-1j * liouvillian * param.time_step),
                                                   eigvectors_liouville))

        # exp_mat = la.expm(-1j * liouvillian * param.time_step)
        # spin1_x = sp.spin1_x
        # spin1_z = sp.spin1_z
        # temp1 = np.real(eigvectors_inv_liouville)
        # temp2 = np.real(exp_mat)
        # temp3 = np.real(eigvectors_liouville)
        start = time.time()
        for bench in range(0, int(1E4)):
            exp_test = la.expm(-1j * liouvillian * param.time_step)
            #test = np.matmul(spin1_x, spin1_z)
            #test = np.matmul(eigvectors_inv[count], np.matmul(spin3_s1_z, eigvectors[count]))
            #propagator[count, :] = np.matmul(eigvectors_inv_liouville, np.matmul(test, eigvectors_liouville))
            #test = np.matmul(temp1, np.matmul(temp2, temp3))
            # eigvectors_liouville = np.kron(eigvectors[count], eigvectors[count])
            #hamiltonian_liouville = kron_rmat_eye(total_hamiltonian, 2 ** num_spins)
        print('la.expm', (time.time() - start))

    return propagator


def basis_transform(eigvectors, eigvectors_inv, matrices):
    """ Transforms basis of inputted list of numpy arrays.
    """

    matrices_out = np.copy(matrices)
    for count in range(0, len(matrices)):
        matrices_out[count] = np.matmul(eigvectors_inv, np.matmul(matrices[count], eigvectors))

    return matrices_out


def density_mat_thermal(hamiltonian_ideal):
    """ Calculate initial thermal density matrix from Boltzmann factors.
    """

    # Calculate energies of idealised Hamiltonian
    energies_ideal = np.diag(hamiltonian_ideal)

    # Calculate initial Zeeman basis density matrix from Boltzmann factors
    boltzmann_factors = np.zeros(len(hamiltonian_ideal))
    for count in range(0, hamiltonian_ideal.shape[0]):
        boltzmann_factors[count] = np.exp(-(sc.Planck * energies_ideal[count]) / (sc.Boltzmann * param.temperature))
    density_mat = (1/np.sum(boltzmann_factors)) * np.diagflat(boltzmann_factors)

    return density_mat


def anisotropy_coefficients(angles):
    """ Calculate time independent coefficients for electron g-anisotropy.
    """

    gx = param.electron_frequency * param.gtensor[0]
    gy = param.electron_frequency * param.gtensor[1]
    gz = param.electron_frequency * param.gtensor[2]

    ca = np.cos(angles[0])
    cb = np.cos(angles[1])
    cg = np.cos(angles[2])
    sa = np.sin(angles[0])
    sb = np.sin(angles[1])
    sg = np.sin(angles[2])

    r11 = ca * cb * cg - sa * sg
    r12 = sa * cb * cg + ca * sg
    r13 = -sb * cg
    r21 = -ca * cb * sg - sa * cg
    r22 = -sa * cb * sg + ca * cg
    r23 = sb * sg
    r31 = ca * sb
    r32 = sa * sb
    r33 = cb

    c0 = 1 / 3 * (gx + gy + gz)
    c1 = 2 * np.sqrt(2) / 3 * (gx * r11 * r31 + gy * r12 * r32 + gz * r13 * r33)
    c2 = 2 * np.sqrt(2) / 3 * (gx * r21 * r31 + gy * r22 * r32 + gz * r23 * r33)
    c3 = 1 / 3 * (gx * (r11 ** 2 - r21 ** 2) + gy * (r12 ** 2 - r22 ** 2) + gz * (r13 ** 2 - r23 ** 2))
    c4 = 2 / 3 * (gx * r11 * r21 + gy * r22 * r12 + gz * r13 * r23)

    return c0, c1, c2, c3, c4


def anisotropy(c0, c1, c2, c3, c4, time):
    """ Calculate time dependent electron g-anisotropy.
    """

    g_anisotropy = c0 + c1 * np.cos(2 * np.pi * param.freq_rotor * time) + \
                   c2 * np.sin(2 * np.pi * param.freq_rotor * time) + \
                   c3 * np.cos(2 * np.pi * param.freq_rotor * time * 2) + \
                   c4 * np.sin(2 * np.pi * param.freq_rotor * time * 2)

    return g_anisotropy


def hyperfine(hyperfine_angles, time):
    """ Calculate time dependent hyperfine, neglecting zy coupling.
    """

    hyperfine_zz = param.hyperfine_coupling * (-0.5 * (np.sin(hyperfine_angles[1]) ** 2) *
                                               np.cos(2 * (2 * np.pi * param.freq_rotor * time + hyperfine_angles[2])) +
                                               np.sqrt(2) * np.sin(hyperfine_angles[1]) *
                                               np.cos(hyperfine_angles[1]) *
                                               np.cos(2 * np.pi * param.freq_rotor * time + hyperfine_angles[2]))

    hyperfine_zx = param.hyperfine_coupling * (-0.5 * np.sin(hyperfine_angles[1]) * np.cos(hyperfine_angles[1]) *
                                               np.cos(2 * np.pi * param.freq_rotor * time + hyperfine_angles[2]) -
                                               (np.sqrt(2) / 4) * (np.sin(hyperfine_angles[1]) ** 2) *
                                               np.cos(2 * (2 * np.pi * param.freq_rotor * time + hyperfine_angles[2])) +
                                               (np.sqrt(2) / 4) * (3 * (np.cos(hyperfine_angles[1]) ** 2) - 1))

    return hyperfine_zz, hyperfine_zx


def kron_rmat_eye(mat, val):
    """ Calculates np.kron(mat, np.eye(val)) of square matrix mat
    """

    n = mat.shape[0]
    out = np.zeros((n, val, n, val), dtype=mat.dtype)
    for i in range(0, val):
        out[:, i, :, i] = mat
    out.shape = (n * val, n * val)

    return out


def kron_eye_rmat(val, mat):
    """ Calculates np.kron(np.eye(val), mat) of square matrix mat
    """

    n = mat.shape[0]
    out = np.zeros((val, n, val, n), dtype=mat.dtype)
    for i in range(0, val):
        out[i, :, i, :] = mat
    out.shape = (n * val, n * val)

    return out
