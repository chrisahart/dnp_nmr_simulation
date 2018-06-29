import parameters as param
import numpy as np
from scipy import linalg as la
from scipy import constants as sc
import spin_matrices as sp
import functions as fn
import time


def dynamics(microwave_amplitude):
    """ Solid effect dynamics for an e-n system.
    """

    # Pre-allocate arrays
    propagator_strobe = np.eye(4 ** param.num_spins)

    # Calculate thermal density matrix from idealised Hamiltonian
    hamiltonian_ideal = param.electron_frequency * sp.spin2_s_z + \
                        param.freq_nuclear_1 * sp.spin2_i_z
    density_mat = fn.density_mat_thermal(hamiltonian_ideal)

    # Construct intrinsic Hilbert space Hamiltonian
    hamiltonian = calculate_hamiltonian()

    # Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
    eigvals, eigvectors = np.linalg.eig(hamiltonian)
    energies = np.real(eigvals)
    eigvectors_inv = np.linalg.inv(eigvectors)

    # Calculate microwave Hamiltonian
    microwave_hamiltonian_init = microwave_amplitude * sp.spin2_s_x

    # Calculate Liouville space propagator with relaxation
    propagator = fn.liouville_propagator(energies, eigvectors, eigvectors_inv, microwave_hamiltonian_init,
                                         relaxation_mat)

    # Calculate stroboscopic propagator (product of all operators within rotor period)
    for count in range(0, int(param.time_step_num)):
        propagator_strobe = np.matmul(propagator_strobe, propagator[count, :])

    # Propagate density matrix, calculating polarisations
    pol_i_z, pol_s_z = calculate_polarisation(density_mat, propagator_strobe)

    return pol_i_z, pol_s_z


def calculate_hamiltonian():

    # Pre-allocate hamiltonian
    hamiltonian = np.zeros((int(param.time_step_num), 2 ** param.num_spins, 2 ** param.num_spins,), dtype=np.complex)

    # Calculate time independent electron g-anisotropy coefficients
    c0, c1, c2, c3, c4 = fn.anisotropy_coefficients(param.orientation_tempol)

    for count in range(0, int(param.time_step_num)):

        # Calculate time dependent hyperfine
        hyperfine_zz, hyperfine_zx = fn.hyperfine(param.hyperfine_angles_1, count * param.time_step)
        hyperfine_total = hyperfine_zz * np.matmul(sp.spin2_i_z, sp.spin2_s_z) + \
                          hyperfine_zx * (np.matmul(sp.spin2_i_p, sp.spin2_s_z) + np.matmul(sp.spin2_i_m, sp.spin2_s_z))

        # Calculate time dependent electron g-anisotropy
        ganisotropy = fn.anisotropy(c0, c1, c2, c3, c4, count * param.time_step)

        # Calculate time dependent hamiltonian
        hamiltonian[count, :] = (ganisotropy - param.microwave_frequency) * sp.spin2_s_z + \
                                param.freq_nuclear_1 * sp.spin2_i_z + \
                                hyperfine_total

    return hamiltonian


def calculate_polarisation(density_mat, propagator_strobe):

    pol_i_z = np.zeros(param.num_timesteps_prop)
    pol_s_z = np.zeros(param.num_timesteps_prop)

    for count in range(0, param.num_timesteps_prop):

        # Calculate electronic and nuclear polarisation
        pol_i_z[count] = np.trace(np.matmul(np.real(density_mat), sp.spin2_i_z))
        pol_s_z[count] = np.trace(np.matmul(np.real(density_mat), sp.spin2_s_z))

        # Transform density matrix (2^N x 2^N to 4^N x 1)
        density_mat_liouville = np.reshape(density_mat, [4 ** param.num_spins, 1])

        # Propagate density matrix
        density_mat = np.matmul(propagator_strobe, density_mat_liouville)

        # Transform density matrix (4^N x 1 to 2^N x 2^N)
        density_mat = np.reshape(density_mat, [2 ** param.num_spins, 2 ** param.num_spins])

    return pol_i_z, pol_s_z


def relaxation_mat(eigvectors, eigvectors_inv, gnp, gnm, gep, gem):
    """ Calculate time dependent Liouville space relaxation matrix.
    Replace np.kron(A, np.eye(2 ** param.num_spins))) with fn.kron_a_n(A, 2 ** param.num_spins)
    """

    # Transform spin matrices into time dependent Hilbert space basis
    spin2_s_x_t, spin2_s_y_t, spin2_s_z_t, spin2_s_p_t, spin2_s_m_t, \
        spin2_i_x_t, spin2_i_y_t, spin2_i_z_t, spin2_i_p_t, spin2_i_m_t = \
            fn.basis_transform(eigvectors, eigvectors_inv, sp.spin2_all)

    # fn.kron_a_n(spin2_i_z_t, 2 ** param.num_spins)

    # Transform spin matrices into time dependent Liouville space basis
    spin2_s_z_tl = np.kron(spin2_s_z_t, np.transpose(spin2_s_z_t))
    spin2_i_z_tl = np.kron(spin2_i_z_t, np.transpose(spin2_i_z_t))
    spin2_i_p_tl = np.kron(spin2_i_p_t, np.transpose(spin2_i_m_t)) - 0.5 * np.eye(4 ** param.num_spins) + 0.5 * (
        np.kron(spin2_i_z_t, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins),
                                                                     np.transpose(spin2_i_z_t)))
    spin2_i_m_tl = np.kron(spin2_i_m_t, np.transpose(spin2_i_p_t)) - 0.5 * np.eye(4 ** param.num_spins) - 0.5 * (
        np.kron(spin2_i_z_t, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins),
                                                                     np.transpose(spin2_i_z_t)))
    spin2_s_p_tl = np.kron(spin2_s_p_t, np.transpose(spin2_s_m_t)) - 0.5 * np.eye(4 ** param.num_spins) + 0.5 * (
        np.kron(spin2_s_z_t, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins),
                                                                     np.transpose(spin2_s_z_t)))
    spin2_s_m_tl = np.kron(spin2_s_m_t, np.transpose(spin2_s_p_t)) - 0.5 * np.eye(4 ** param.num_spins) - 0.5 * (
        np.kron(spin2_s_z_t, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins),
                                                                     np.transpose(spin2_s_z_t)))

    # Calculate relaxation matrices
    relax_t2_elec = (1 / param.t2_elec) * (spin2_s_z_tl - 0.25 * np.eye(4 ** param.num_spins))
    relax_t2_nuc = (1 / param.t2_nuc) * (spin2_i_z_tl - 0.25 * np.eye(4 ** param.num_spins))
    relax_t1 = gep * spin2_s_p_tl + gem * spin2_s_m_tl + gnp * spin2_i_p_tl + gnm * spin2_i_m_tl
    relax_mat = relax_t2_elec + relax_t2_nuc + relax_t1

    return relax_mat
