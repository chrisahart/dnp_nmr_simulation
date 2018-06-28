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

    start = time.time()

    # Pre-allocate arrays for Hamiltonian and propagator
    hamiltonian = np.zeros((int(param.time_step_num), 2 ** param.num_spins, 2 ** param.num_spins,), dtype=np.complex)
    propagator = np.zeros((int(param.time_step_num), 4 ** param.num_spins, 4 ** param.num_spins), dtype=np.complex)
    propagator_strobe = np.eye(4 ** param.num_spins)

    # Pre-allocate arrays for polarisation
    pol_i_z = np.zeros(param.num_timesteps_prop)
    pol_s_z = np.zeros(param.num_timesteps_prop)

    # Calculate variables for Liouville space relaxation
    p_e = np.tanh(0.5 * param.electron_frequency * (sc.Planck / (sc.Boltzmann * param.temperature)))
    p_n = np.tanh(0.5 * param.freq_nuclear_1 * (sc.Planck / (sc.Boltzmann * param.temperature)))
    gnp = 0.5 * (1 - p_n) * (1 / (1 * param.t1_nuc))
    gnm = 0.5 * (1 + p_n) * (1 / (1 * param.t1_nuc))
    gep = 0.5 * (1 - p_e) * (1 / (1 * param.t1_elec))
    gem = 0.5 * (1 + p_e) * (1 / (1 * param.t1_elec))

    # Calculate time independent electron g-anisotropy coefficients
    c0, c1, c2, c3, c4 = fn.anisotropy_coefficients(param.orientation_tempol)

    # Calculate thermal density matrix from idealised Hamiltonian
    hamiltonian_ideal = param.electron_frequency * sp.spin2_s_z + \
                        param.freq_nuclear_1 * sp.spin2_i_z
    density_mat = fn.density_mat_thermal(hamiltonian_ideal)

    print('Density matrix', time.time() - start)

    # Construct intrinsic Hamiltonian (Hilbert space)
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

    print('Hamiltonian', time.time() - start)

    # Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
    eigvals, eigvectors = np.linalg.eig(hamiltonian)
    energies = np.real(eigvals)
    eigvectors_inv = np.linalg.inv(eigvectors)

    # Calculate microwave Hamiltonian
    microwave_hamiltonian_init = microwave_amplitude * sp.spin2_s_x

    print('Eigenvalues', time.time() - start)

    # Calculate Liouville space propagator with relaxation todo parallelise? each loop takes 2ms with 10000 loops
    for count in range(0, int(param.time_step_num)):

        loop = time.time()
        print('number of loops', int(param.time_step_num))

        # Transform microwave Hamiltonian into time dependent basis
        microwave_hamiltonian = np.matmul(eigvectors_inv[count], np.matmul(microwave_hamiltonian_init,
                                                                           eigvectors[count]))

        # Calculate total Hamiltonian
        total_hamiltonian = np.diagflat(energies[count]) + microwave_hamiltonian

        # Transform Hilbert space Hamiltonian into Liouville space
        hamiltonian_liouville = np.kron(total_hamiltonian, np.eye(4)) - np.kron(np.eye(4),
                                                                                np.transpose(total_hamiltonian))

        # Calculate time dependent Liouville space relaxation matrix
        relax_mat = relaxation_mat(eigvectors[count], eigvectors_inv[count],
                                   gnp, gnm, gep, gem)

        # Calculate Liouville space eigenvectors
        eigvectors_liouville = np.kron(eigvectors[count], eigvectors[count])
        eigvectors_inv_liouville = np.kron(eigvectors_inv[count], eigvectors_inv[count])

        # Calculate Liouville space propagator
        liouvillian = hamiltonian_liouville + 1j * relax_mat
        propagator[count, :] = np.matmul(eigvectors_inv_liouville,
                                         np.matmul(la.expm(-1j * liouvillian * param.time_step), eigvectors_liouville))

        print('loop', time.time() - loop)

    print('Calculate Liouville space propagator', time.time() - start)

    # Calculate stroboscopic propagator (product of all operators within rotor period)
    for count in range(0, int(param.time_step_num)):
        propagator_strobe = np.matmul(propagator_strobe, propagator[count, :])

    print('Calculate stroboscopic propagator', time.time() - start)

    # Propagate density matrix, calculating polarisations
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

    print('Propagate density matrix', time.time() - start)

    return pol_i_z, pol_s_z


def relaxation_mat(eigvectors, eigvectors_inv, gnp, gnm, gep, gem):
    """ Calculate time dependent Liouville space relaxation matrix.
    """

    # Transform spin matrices into time dependent Hilbert space basis
    spin2_s_x_t, spin2_s_y_t, spin2_s_z_t, spin2_s_p_t, spin2_s_m_t, \
        spin2_i_x_t, spin2_i_y_t, spin2_i_z_t, spin2_i_p_t, spin2_i_m_t = \
            fn.basis_transform(eigvectors, eigvectors_inv, sp.spin2_all)

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
