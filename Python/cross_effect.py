import parameters as param
import numpy as np
from scipy import linalg as la
from scipy import constants as sc
import spin_matrices as sp
import functions as fn
import time
import cross_effect_plotting
from shutil import copyfile
import os
import matplotlib.pyplot as plt
import f2py_cross_effect as fortran


def main():
    """ Loop over dynamics function as required, writing output to disk and calling plotting function.
    """

    # Pre-allocate arrays
    pol_nuc = np.zeros((param.microwave_amplitude.size, int(param.num_timesteps_prop)))
    pol_elec1 = np.zeros((param.microwave_amplitude.size, int(param.num_timesteps_prop)))
    pol_elec2 = np.zeros((param.microwave_amplitude.size, int(param.num_timesteps_prop)))

    # Start timer
    start = time.time()

    # Calculate system dynamics, looping over microwave amplitude array
    for count in range(0, param.microwave_amplitude.size):

        pol_nuc[count], pol_elec1[count], pol_elec2[count], pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot, energies = \
            dynamics(param.microwave_amplitude[count])

        print('{}{:d}{}{:d}{}{:.2f}{}'.format('Finished loop ', (count + 1), ' of ', param.microwave_amplitude.size,
                                              ', total elapsed time ', (time.time() - start), ' s.'))

    # # # Dynamically assign and create output directory
    # directory = '{}{:.2f}'.format('out/cross_effect/mw_', param.microwave_amplitude[-1] / 1E6)
    # if not os.path.exists(directory):
    #     os.makedirs(directory)
    #
    # # Save data and copy parameters file
    # np.savetxt('{}{}'.format(directory, '/pol_nuc.csv'), pol_nuc, fmt='%.8f', newline='\n')
    # np.savetxt('{}{}'.format(directory, '/pol_elec1.csv'), pol_elec1, fmt='%.8f', newline='\n')
    # np.savetxt('{}{}'.format(directory, '/pol_elec2.csv'), pol_elec2, fmt='%.8f', newline='\n')
    #
    # # High precision required for sub rotor dynamics (negligible change in file size)
    # np.savetxt('{}{}'.format(directory, '/pol_i_z_rot.csv'), pol_i_z_rot, fmt='%.12f', newline='\n')
    # np.savetxt('{}{}'.format(directory, '/pol_s1_z_rot.csv'), pol_s1_z_rot, fmt='%.12f', newline='\n')
    # np.savetxt('{}{}'.format(directory, '/pol_s2_z_rot.csv'), pol_s2_z_rot, fmt='%.12f', newline='\n')
    # np.savetxt('{}{}'.format(directory, '/energies.csv'), energies, fmt='%.0f', newline='\n')
    #
    # # Create a copy of parameters file
    # copyfile("parameters.py", '{}{}'.format(directory, '/parameters.py'))
    #
    # # Call plotting function
    # cross_effect_plotting.plot_all(directory)


def dynamics(microwave_amplitude):
    """ Calculate dynamics of e-n system, calling required functions.
    """

    # Pre-allocate arrays
    propagator_strobe = np.eye(4 ** 3)

    # Pre-compute spin matrices to avoid re-evaluating Kronecker products
    spin3_s1_x = sp.spin3_s1_x
    spin3_s1_z = sp.spin3_s1_z
    spin3_s1_p = sp.spin3_s1_p
    spin3_s1_m = sp.spin3_s1_m

    # 8x8 Matrix operators for S2 operator
    spin3_s2_x = sp.spin3_s2_x
    spin3_s2_z = sp.spin3_s2_z
    spin3_s2_p = sp.spin3_s2_p
    spin3_s2_m = sp.spin3_s2_m

    # 8x8 Matrix operators for I operator
    spin3_i_x = sp.spin3_i_x
    spin3_i_z = sp.spin3_i_z
    spin3_all = sp.spin3_all

    # Calculate thermal density matrix from idealised Hamiltonian
    hamiltonian_ideal = param.electron_frequency * spin3_s1_z + \
                        param.electron_frequency * spin3_s2_z + \
                        param.freq_nuclear_1 * spin3_i_z
    density_mat = fn.density_mat_thermal(hamiltonian_ideal)

    # Construct intrinsic Hilbert space Hamiltonian
    hamiltonian = calculate_hamiltonian(spin3_s1_z, spin3_s2_z, spin3_i_x, spin3_i_z,
                                        spin3_s1_m, spin3_s1_p, spin3_s2_m, spin3_s2_p)

    # Calculate thermal density matrix and intrinsic Hilbert space Hamiltonian using F2PY
    # hamiltonian, density_mat = fortran.cross_effect_dynamics.calculate_hamiltonian(
    #     param.time_step_num, param.time_step, param.freq_rotor, param.gtensor, param.temperature,
    #     param.hyperfine_coupling, param.hyperfine_angles_1, param.orientation_ce_1, param.orientation_ce_2,
    #     param.electron_frequency, param.microwave_frequency, param.freq_nuclear_1, 8, 64)

    # Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
    eigvals, eigvectors = np.linalg.eig(hamiltonian)
    energies = np.copy(eigvals)

    # Determine indices for sorting eigenvalues and eigenvectors into descending order at zero time
    indices = np.argsort(-eigvals[0])

    # Sort eigenvalues and eigenvectors at zero time
    eigvectors_temp = np.copy(eigvectors)
    for count2 in range(0, param.time_step_num):
        for count3 in range(0, 2 ** 3):
            energies[count2, count3] = eigvals[count2, indices[count3]]
            eigvectors[count2, :, count3] = eigvectors_temp[count2, :, indices[count3]]

    # Calculate inverse eigenvectors
    eigvectors_inv = np.linalg.inv(eigvectors)

    # Calculate microwave Hamiltonian
    microwave_hamiltonian_init = microwave_amplitude * (spin3_s1_x + spin3_s2_x)

    # Calculate Liouville space propagator with relaxation
    propagator = fn.liouville_propagator(3, energies, eigvectors, eigvectors_inv,
                                         microwave_hamiltonian_init, calculate_relaxation_mat, spin3_all)

    # Calculate Liouville space propagator with relaxation using F2PY
    # propagator = fortran.cross_effect_dynamics.liouville_propagator(
    #     param.time_step_num, param.time_step, param.electron_frequency, param.freq_nuclear_1,
    #     param.microwave_amplitude, param.t1_nuc, param.t1_elec, param.t2_nuc, param.t2_elec, param.temperature,
    #     8, 64, eigvectors, eigvectors_inv, energies)

    # Propagate density matrix for single rotor period, calculating polarisations
    pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot = calculate_polarisation_sub_rotor(density_mat, propagator,
                                                                                spin3_s1_z, spin3_s2_z, spin3_i_z)

    # Propagate density matrix for single rotor period, calculating polarisations using F2PY
    # pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot = fortran.cross_effect_dynamics.calculate_polarisation_sub_rotor(
    #     param.time_step_num, density_mat, propagator, 8, 64)

    # Calculate stroboscopic propagator (product of all operators within rotor period)
    for count in range(0, int(param.time_step_num)):
        propagator_strobe = np.matmul(propagator_strobe, propagator[count, :])

    # Propagate density matrix stroboscopically, calculating polarisations
    pol_i_z, pol_s1_z, pol_s2_z = calculate_polarisation_rotor(density_mat, propagator_strobe,
                                                               spin3_s1_z, spin3_s2_z, spin3_i_z)

    # Propagate density matrix stroboscopically, calculating polarisations using F2PY
    # pol_i_z, pol_s1_z, pol_s2_z = fortran.cross_effect_dynamics.calculate_polarisation_rotor(
    #     param.time_step_num, param.num_timesteps_prop, density_mat, propagator, 8, 64)

    return pol_i_z, pol_s1_z, pol_s2_z, pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot, energies


def calculate_hamiltonian(spin3_s1_z, spin3_s2_z, spin3_i_x, spin3_i_z,
                          spin3_s1_m, spin3_s1_p, spin3_s2_m, spin3_s2_p):
    """ Calculate Hamiltonian of e-n system.
    """

    # Pre-allocate hamiltonian
    hamiltonian = np.zeros((int(param.time_step_num), 2 ** 3, 2 ** 3,))

    # Calculate time independent electron g-anisotropy coefficients for electron 1 and 2
    c0_1, c1_1, c2_1, c3_1, c4_1 = fn.anisotropy_coefficients(param.orientation_ce_1)
    c0_2, c1_2, c2_2, c3_2, c4_2 = fn.anisotropy_coefficients(param.orientation_ce_2)

    for count in range(0, int(param.time_step_num)):

        # Calculate time dependent hyperfine between electron 1 and nucleus
        hyperfine_zz, hyperfine_zx = fn.hyperfine(param.hyperfine_angles_1, count * param.time_step)
        hyperfine_total = hyperfine_zz * 2 * np.matmul(spin3_i_z, spin3_s1_z) + \
                          hyperfine_zx * 2 * np.matmul(spin3_i_x, spin3_s1_z)

        # Calculate time dependent electron g-anisotropy for electron 1 and 2
        ganisotropy_1 = fn.anisotropy(c0_1, c1_1, c2_1, c3_1, c4_1, count * param.time_step)
        ganisotropy_2 = fn.anisotropy(c0_2, c1_2, c2_2, c3_2, c4_2, count * param.time_step)

        # Calculate time dependent dipolar between electron 1 and 2
        # S20 = 2 * np.matmul(spin3_s1_z, spin3_s2_z) - \
        #          0.5 * (np.matmul(spin3_s1_p, spin3_s2_m) + np.matmul(spin3_s1_m, spin3_s2_p))
        # b_ee = 0
        # g_ee = 0
        # test = 23e6 * (+0.5 * ((np.sin(b_ee)) ** 2) * np.cos(2 * 2 * np.pi * param.freq_rotor *
        #                                                     count * param.time_step_num + g_ee) -
        #                np.sqrt(2) * np.sin(b_ee) * np.cos(b_ee) * np.cos(2 * np.pi * param.freq_rotor *
        #                                                                count * param.time_step_num + g_ee))
        # dipolar = 0  # test * S20

        # Calculate time dependent hamiltonian
        hamiltonian[count, :] = (ganisotropy_1 - param.microwave_frequency) * spin3_s1_z + \
                                (ganisotropy_2 - param.microwave_frequency) * spin3_s2_z + \
                                param.freq_nuclear_1 * spin3_i_z + \
                                hyperfine_total

    return hamiltonian


def calculate_polarisation_rotor(density_mat, propagator_strobe, spin3_s1_z, spin3_s2_z, spin3_i_z):
    """ Calculate polarisation across multiple rotor periods.
    """

    pol_i_z = np.zeros(param.num_timesteps_prop)
    pol_s1_z = np.zeros(param.num_timesteps_prop)
    pol_s2_z = np.zeros(param.num_timesteps_prop)

    for count in range(0, param.num_timesteps_prop):

        # Calculate electronic and nuclear polarisation
        pol_i_z[count] = np.trace(np.matmul(np.real(density_mat), spin3_i_z))
        pol_s1_z[count] = np.trace(np.matmul(np.real(density_mat), spin3_s1_z))
        pol_s2_z[count] = np.trace(np.matmul(np.real(density_mat), spin3_s2_z))

        # Transform density matrix (2^N x 2^N to 4^N x 1)
        density_mat_liouville = np.reshape(density_mat, [4 ** 3, 1])

        density_mat = np.matmul(propagator_strobe, density_mat_liouville)

        # Transform density matrix (4^N x 1 to 2^N x 2^N)
        density_mat = np.reshape(density_mat, [2 ** 3, 2 ** 3])

    return pol_i_z, pol_s1_z, pol_s2_z


def calculate_polarisation_sub_rotor(density_mat, propagator, spin3_s1_z, spin3_s2_z, spin3_i_z):
    """ Calculate polarisation across single rotor period.
    """

    pol_i_z_rot = np.zeros(param.time_step_num)
    pol_s1_z_rot = np.zeros(param.time_step_num)
    pol_s2_z_rot = np.zeros(param.time_step_num)

    for count in range(0, param.time_step_num):

        # Calculate electronic and nuclear polarisation
        pol_i_z_rot[count] = np.trace(np.matmul(np.real(density_mat), spin3_i_z))
        pol_s1_z_rot[count] = np.trace(np.matmul(np.real(density_mat), spin3_s1_z))
        pol_s2_z_rot[count] = np.trace(np.matmul(np.real(density_mat), spin3_s2_z))

        # Transform density matrix (2^N x 2^N to 4^N x 1)
        density_mat_liouville = np.reshape(density_mat, [4 ** 3, 1])

        # Propagate density matrix
        density_mat = np.matmul(propagator[count, :], density_mat_liouville)

        # Transform density matrix (4^N x 1 to 2^N x 2^N)
        density_mat = np.reshape(density_mat, [2 ** 3, 2 ** 3])

    return pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot


def calculate_relaxation_mat(eigvectors, eigvectors_inv, gnp, gnm, gep, gem, spin3_all):
    """ Calculate time dependent Liouville space relaxation matrix.
    """

    identity_mat = 0.5 * np.eye(4 ** 3)

    # Transform spin matrices into time dependent Hilbert space basis
    spin3_s1_z_t, spin3_s1_x_t, spin3_s1_p_t, spin3_s1_m_t, spin3_s2_z_t, spin3_s2_x_t, spin3_s2_p_t, spin3_s2_m_t, \
        spin3_i_z_t, spin3_i_x_t, spin3_i_p_t, spin3_i_m_t = \
        fn.basis_transform(eigvectors, eigvectors_inv, spin3_all)

    # Transform spin matrices into time dependent Liouville space basis
    spin3_i_p_tl = np.kron(spin3_i_p_t, np.transpose(spin3_i_m_t)) - identity_mat + 0.5 * (
            fn.kron_rmat_eye(spin3_i_z_t, 8) + fn.kron_eye_rmat(8, np.transpose(spin3_i_z_t)))

    spin3_i_m_tl = np.kron(spin3_i_m_t, np.transpose(spin3_i_p_t)) - identity_mat - 0.5 * (
            fn.kron_rmat_eye(spin3_i_z_t, 8) + fn.kron_eye_rmat(8, np.transpose(spin3_i_z_t)))

    spin3_s1_p_tl = np.kron(spin3_s1_p_t, np.transpose(spin3_s1_m_t)) - identity_mat + 0.5 * (
            fn.kron_rmat_eye(spin3_s1_z_t, 8) + fn.kron_eye_rmat(8, np.transpose(spin3_s1_z_t)))

    spin3_s1_m_tl = np.kron(spin3_s1_m_t, np.transpose(spin3_s1_p_t)) - identity_mat - 0.5 * (
            fn.kron_rmat_eye(spin3_s1_z_t, 8) + fn.kron_eye_rmat(8, np.transpose(spin3_s1_z_t)))

    spin3_s2_p_tl = np.kron(spin3_s2_p_t, np.transpose(spin3_s2_m_t)) - identity_mat + 0.5 * (
            fn.kron_rmat_eye(spin3_s2_z_t, 8) + fn.kron_eye_rmat(8, np.transpose(spin3_s2_z_t)))

    spin3_s2_m_tl = np.kron(spin3_s2_m_t, np.transpose(spin3_s2_p_t)) - identity_mat - 0.5 * (
            fn.kron_rmat_eye(spin3_s2_z_t, 8) + fn.kron_eye_rmat(8, np.transpose(spin3_s2_z_t)))

    # Calculate relaxation matrices
    relax_t2_s1 = (2 / param.t2_elec) * (np.kron(spin3_s1_z_t, np.transpose(spin3_s1_z_t)) - 0.5 * identity_mat)
    relax_t2_s2 = (2 / param.t2_elec) * (np.kron(spin3_s2_z_t, np.transpose(spin3_s2_z_t)) - 0.5 * identity_mat)
    relax_t2_i1 = (2 / param.t2_nuc) * (np.kron(spin3_i_z_t, np.transpose(spin3_i_z_t)) - 0.5 * identity_mat)
    relax_t1 = gep * spin3_s1_p_tl + gem * spin3_s1_m_tl + \
               gep * spin3_s2_p_tl + gem * spin3_s2_m_tl + \
               gnp * spin3_i_p_tl + gnm * spin3_i_m_tl
    relax_mat = relax_t2_s1 + relax_t2_s2 + relax_t2_i1 + relax_t1

    return relax_mat


def calculate_relaxation_mat_mance(eigvectors, eigvectors_inv, gnp, gnm, gep, gem, spin3_all):
    """ Calculate time dependent Liouville space relaxation matrix from Mance paper.
    """

    # Transform spin matrices into time dependent Hilbert space basis
    spin3_s1_z_t, spin3_s1_x_t, spin3_s1_p_t, spin3_s1_m_t, spin3_s2_z_t, spin3_s2_x_t, spin3_s2_p_t, spin3_s2_m_t, \
        spin3_i_z_t, spin3_i_x_t, spin3_i_p_t, spin3_i_m_t = \
        fn.basis_transform(eigvectors, eigvectors_inv, spin3_all)

    # Calculate Boltzmann factors and pre-allocate variables
    boltzmann_elec = param.electron_frequency * (sc.Planck / (sc.Boltzmann * param.temperature))
    boltzmann_nuc = param.freq_nuclear_1 * (sc.Planck / (sc.Boltzmann * param.temperature))
    relax_values_t1 = np.zeros((2 ** 3, 2 ** 3))
    relax_mat_t1 = np.zeros((4 ** 3, 4 ** 3))

    # Calculate T1 relaxation matrix
    for a in range(0, 2 ** 3):
        for b in range(0, 2 ** 3):

            # Calculate Boltzmann factors
            diff_elec = spin3_s1_z_t[b, b] - spin3_s1_z_t[a, a] + spin3_s2_z_t[b, b] - spin3_s2_z_t[a, a]
            diff_nuc = spin3_i_z_t[b, b] - spin3_i_z_t[a, a]
            boltzmann_factor = np.exp(-diff_elec * boltzmann_elec / 2 - diff_nuc * boltzmann_nuc / 2) / \
                               (np.exp(diff_elec * boltzmann_elec / 2 + diff_nuc * boltzmann_nuc / 2) +
                                np.exp(-diff_elec * boltzmann_elec / 2 - diff_nuc * boltzmann_nuc / 2))

            # Calculate initial relaxation matrix
            relax_values_t1[b, a] = (1 / param.t1_elec) * (spin3_s1_x_t[b, a] * spin3_s1_x_t[a, b] +
                                                           spin3_s2_x_t[b, a] * spin3_s2_x_t[a, b] +
                                                           spin3_s1_z_t[b, a] * spin3_s1_z_t[a, b] +
                                                           spin3_s2_z_t[b, a] * spin3_s2_z_t[a, b]) + \
                                    (1 / param.t1_nuc) * (spin3_i_x_t[b, a] * spin3_i_x_t[a, b] +
                                                          spin3_i_z_t[b, a] * spin3_i_z_t[a, b])
            relax_values_t1[b, a] = relax_values_t1[b, a] * boltzmann_factor

    # Set diagonal elements
    np.fill_diagonal(relax_values_t1, 0)
    for a in range(0, 2 ** 3):
        for b in range(0, 2 ** 3):
            if abs(a - b) > 0:
                relax_values_t1[b, b] = relax_values_t1[b, b] - relax_values_t1[a, b]

    # Transform to 4 ** N by 4 ** N matrix
    for a in range(0, 2 ** 3):
        for b in range(0, 2 ** 3):
            c = a * (2 ** 3) + a
            d = b * (2 ** 3) + b
            relax_mat_t1[d, c] = relax_values_t1[b, a]

    # Calculate T2 relaxation matrix using circular shift of spin matrices
    relax_values_t2 = np.zeros(2 ** 3)
    for a in range(1, 2 ** 3):
        relax_values_t2[a] = a
        relax_values_t2[a] = (((abs(spin3_s1_z_t[a, a] - spin3_s1_z_t[a - 1, a - 1])) ** 2) +
                              ((abs(spin3_s2_z_t[a, a] - spin3_s2_z_t[a - 1, a - 1])) ** 2)) \
                             * (1 / param.t2_elec) + \
                             ((abs(spin3_i_z_t[a, a] - spin3_i_z_t[a - 1, a - 1])) ** 2) * (1 / param.t2_nuc)

    # Transform to 4 ** N by 4 ** N matrix
    relax_mat_t2 = np.diagflat(np.concatenate((relax_values_t2, np.roll(relax_values_t2, 1),
                                               np.roll(relax_values_t2, 2), np.roll(relax_values_t2, 3),
                                               np.roll(relax_values_t2, 4), np.roll(relax_values_t2, 5),
                                               np.roll(relax_values_t2, 6), np.roll(relax_values_t2, 7))))

    # Relaxation matrix as sum of t1 and t2 matrices
    relax_mat = -1 * relax_mat_t2 + relax_mat_t1

    return relax_mat


if __name__ == "__main__":
    main()
    print('End of program.')
    plt.show()
