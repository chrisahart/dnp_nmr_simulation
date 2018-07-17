import parameters as param
import numpy as np
from scipy import linalg as la
from scipy import constants as sc
import spin_matrices as sp
import functions as fn
import time
import solid_effect_plotting
from shutil import copyfile
import os
import matplotlib.pyplot as plt
import f2py_dynamics as fortran


def main():
    """ Loop over dynamics function as required, writing output to disk and calling plotting function.
    """

    # Pre-allocate arrays
    pol_nuc = np.zeros((param.microwave_amplitude.size, int(param.num_timesteps_prop)))
    pol_elec = np.zeros((param.microwave_amplitude.size, int(param.num_timesteps_prop)))

    # Start timer
    start = time.time()

    # Calculate system dynamics, looping over microwave amplitude array
    for count in range(0, param.microwave_amplitude.size):
        pol_nuc[count], pol_elec[count], pol_i_z_rot, pol_s_z_rot, energies = \
            dynamics(param.microwave_amplitude[count])

        print('{}{:d}{}{:d}{}{:.2f}{}'.format('Finished loop ', (count + 1), ' of ', param.microwave_amplitude.size,
                                              ', total elapsed time ', (time.time() - start), ' s.'))

    # Dynamically assign and create output directory
    directory = '{}{:.2f}'.format('out/solid_effect/mw_', param.microwave_amplitude[-1] / 1E6)
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Save data and copy parameters file
    np.savetxt('{}{}'.format(directory, '/pol_nuc.csv'), pol_nuc, fmt='%.8f', newline='\n')
    np.savetxt('{}{}'.format(directory, '/pol_elec.csv'), pol_elec, fmt='%.8f', newline='\n')

    # High precision required for sub rotor dynamics (negligible change in file size)
    np.savetxt('{}{}'.format(directory, '/pol_i_z_rot.csv'), pol_i_z_rot, fmt='%.12f', newline='\n')
    np.savetxt('{}{}'.format(directory, '/pol_s_z_rot.csv'), pol_s_z_rot, fmt='%.12f', newline='\n')
    np.savetxt('{}{}'.format(directory, '/energies.csv'), energies, fmt='%.0f', newline='\n')

    # Create a copy of parameters file
    copyfile("parameters.py", '{}{}'.format(directory, '/parameters.py'))

    # Call plotting function
    solid_effect_plotting.plot_all(directory)


def dynamics(microwave_amplitude):
    """ Calculate dynamics of e-n system.
        F2PY Fortran code can be run via single main() call, or per function for comparison purposes.
        Python code is left with some optimisation omitted for clarity, as performance gain is negligible.
    """

    # Pre-compute spin matrices to avoid re-evaluation Kronecker products
    spin2_s_x = sp.spin2_s_x
    spin2_s_z = sp.spin2_s_z

    # Pre-compute spin matrices to avoid re-evaluation Kronecker products
    spin2_i_x = sp.spin2_i_x
    spin2_i_z = sp.spin2_i_z
    spin2_all = sp.spin2_all

    # Calculate thermal density matrix from idealised Hamiltonian
    hamiltonian_ideal = param.electron_frequency * spin2_s_z + \
                        param.freq_nuclear_1 * spin2_i_z
    density_mat = fn.density_mat_thermal(hamiltonian_ideal)

    # Construct intrinsic Hilbert space Hamiltonian
    hamiltonian = calculate_hamiltonian(spin2_s_z, spin2_i_x, spin2_i_z)
    # energies, eigvectors, eigvectors_inv, density_mat = fortran.f2py_dynamics.calculate_hamiltonian(
    # param.time_step_num,
    #                                                                                 param.time_step, param.freq_rotor,
    #                                                          param.gtensor,
    #                                                          param.hyperfine_coupling, param.hyperfine_angles_1,
    #                                                          param.orientation_se, param.electron_frequency,
    #                                                          param.microwave_frequency, param.freq_nuclear_1)

    # Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
    eigvals, eigvectors = np.linalg.eig(hamiltonian)
    energies = np.copy(eigvals)

    # Determine indices for sorting eigenvalues and eigenvectors into descending order at zero time
    indices = np.argsort(-eigvals[0])

    # Sort eigenvalues and eigenvectors at zero time
    eigvectors_temp = np.copy(eigvectors)
    for count2 in range(0, param.time_step_num):
        for count3 in range(0, 4):
            energies[count2, count3] = eigvals[count2, indices[count3]]
            eigvectors[count2, :, count3] = eigvectors_temp[count2, :, indices[count3]]

    # Calculate inverse eigenvectors
    eigvectors_inv = np.linalg.inv(eigvectors)

    # Calculate microwave Hamiltonian
    microwave_hamiltonian_init = microwave_amplitude * spin2_s_x

    # Calculate Liouville space propagator with relaxation
    propagator = fn.liouville_propagator(2, energies, eigvectors, eigvectors_inv,
                                         microwave_hamiltonian_init, calculate_relaxation_mat_mance, spin2_all)
    # propagator = fortran.f2py_dynamics.liouville_propagator(param.time_step_num, param.time_step,
    #                                                         param.electron_frequency, param.freq_nuclear_1,
    #                                                         param.microwave_amplitude, param.t1_nuc,
    #                                                         param.t1_elec, param.t2_nuc, param.t2_elec,
    #                                                         param.temperature, eigvectors,
    #                                                         eigvectors_inv, energies)

    # Propagate density matrix for single rotor period, calculating polarisations
    pol_i_z_rot, pol_s_z_rot = calculate_polarisation_sub_rotor(density_mat, propagator, spin2_s_z, spin2_i_z)
    # pol_i_z_rot, pol_s_z_rot = fortran.f2py_dynamics.calculate_polarisation_sub_rotor(param.time_step_num,
    #                                                                       density_mat, propagator)

    # Calculate stroboscopic propagator (product of all operators within rotor period)
    propagator_strobe = np.eye(4 ** 2)
    for count in range(0, int(param.time_step_num)):
        propagator_strobe = np.matmul(propagator_strobe, propagator[count, :])

    # Propagate density matrix stroboscopically, calculating polarisations
    pol_i_z, pol_s_z = calculate_polarisation_rotor(density_mat, propagator_strobe, spin2_s_z, spin2_i_z)
    # pol_i_z, pol_s_z = fortran.f2py_dynamics.calculate_polarisation_rotor(param.time_step_num,
    #                                                                       param.num_timesteps_prop,
    #                                                                               density_mat, propagator)

    # Call F2PY code to calculate dynamics
    # start = time.time()
    # pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies = fortran.f2py_dynamics.main(param.time_step_num,
    #                                                                                   param.time_step,
    #             param.freq_rotor, param.gtensor, param.hyperfine_coupling, param.hyperfine_angles_1,
    #             param.orientation_se, param.electron_frequency, param.microwave_frequency, param.freq_nuclear_1,
    #                             microwave_amplitude, param.t1_nuc, param.t1_elec, param.t2_nuc, param.t2_elec,
    #             param.temperature, param.num_timesteps_prop)
    # end = time.time() - start
    # print('main() time taken', end)

    return pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies


def calculate_hamiltonian(spin2_s_z, spin2_i_x, spin2_i_z):
    """ Calculate Hamiltonian of e-n system. Use of functions makes negligible difference to benchmark.
    """

    # Pre-allocate hamiltonian
    hamiltonian = np.zeros((int(param.time_step_num), 2 ** 2, 2 ** 2))

    # Calculate time independent electron g-anisotropy coefficients
    c0, c1, c2, c3, c4 = fn.anisotropy_coefficients(param.orientation_se)

    for count in range(0, int(param.time_step_num)):

        # Calculate time dependent hyperfine
        hyperfine_zz, hyperfine_zx = fn.hyperfine(param.hyperfine_angles_1, count * param.time_step)
        hyperfine_total = hyperfine_zz * 2 * np.matmul(spin2_i_z, spin2_s_z) + \
                          hyperfine_zx * 2 * np.matmul(spin2_i_x, spin2_s_z)

        # Calculate time dependent electron g-anisotropy
        ganisotropy = fn.anisotropy(c0, c1, c2, c3, c4, count * param.time_step)

        # Calculate time dependent hamiltonian
        hamiltonian[count, :] = (ganisotropy - param.microwave_frequency) * spin2_s_z + \
                                param.freq_nuclear_1 * spin2_i_z + \
                                hyperfine_total

    return hamiltonian


def calculate_polarisation_rotor(density_mat, propagator_strobe, spin2_s_z, spin2_i_z):
    """ Calculate polarisation across multiple rotor periods.
    """

    pol_i_z = np.zeros(param.num_timesteps_prop)
    pol_s_z = np.zeros(param.num_timesteps_prop)

    for count in range(0, param.num_timesteps_prop):

        # Calculate electronic and nuclear polarisation
        pol_i_z[count] = np.trace(np.matmul(np.real(density_mat), spin2_i_z))
        pol_s_z[count] = np.trace(np.matmul(np.real(density_mat), spin2_s_z))

        # Transform density matrix (2^N x 2^N to 4^N x 1)
        density_mat_liouville = np.reshape(density_mat, [4 ** 2, 1])

        # Propagate density matrix
        density_mat = np.matmul(propagator_strobe, density_mat_liouville)

        # Transform density matrix (4^N x 1 to 2^N x 2^N)
        density_mat = np.reshape(density_mat, [2 ** 2, 2 ** 2])

    return pol_i_z, pol_s_z


def calculate_polarisation_sub_rotor(density_mat, propagator, spin2_s_z, spin2_i_z):
    """ Calculate polarisation across single rotor period.
    """

    pol_i_z_rot = np.zeros(param.time_step_num)
    pol_s_z_rot = np.zeros(param.time_step_num)

    for count in range(0, param.time_step_num):

        # Calculate electronic and nuclear polarisation
        pol_i_z_rot[count] = np.trace(np.matmul(np.real(density_mat), spin2_i_z))
        pol_s_z_rot[count] = np.trace(np.matmul(np.real(density_mat), spin2_s_z))

        # Transform density matrix (2^N x 2^N to 4^N x 1)
        density_mat_liouville = np.reshape(density_mat, [4 ** 2, 1])

        # Propagate density matrix
        density_mat = np.matmul(propagator[count, :], density_mat_liouville)

        # Transform density matrix (4^N x 1 to 2^N x 2^N)
        density_mat = np.reshape(density_mat, [2 ** 2, 2 ** 2])

    return pol_i_z_rot, pol_s_z_rot


def calculate_relaxation_mat(eigvectors, eigvectors_inv, gnp, gnm, gep, gem, spin2_all):
    """ Calculate time dependent Liouville space relaxation matrix.
    """

    identity_mat = 0.5 * np.eye(4 ** 2)

    # Transform spin matrices into time dependent Hilbert space basis
    spin2_s_z_t, spin2_s_x_t, spin2_s_p_t, spin2_s_m_t, spin2_i_z_t, spin2_i_x_t, spin2_i_p_t, spin2_i_m_t = \
        fn.basis_transform(eigvectors, eigvectors_inv, spin2_all)

    # Transform spin matrices into time dependent Liouville space basis
    spin2_i_p_tl = np.kron(spin2_i_p_t, np.transpose(spin2_i_m_t)) - identity_mat + 0.5 * (
        fn.kron_a_n(spin2_i_z_t, 2 ** 2) +
        fn.kron_n_a(2 ** 2, np.transpose(spin2_i_z_t)))

    spin2_i_m_tl = np.kron(spin2_i_m_t, np.transpose(spin2_i_p_t)) - identity_mat - 0.5 * (
        fn.kron_a_n(spin2_i_z_t, 2 ** 2) +
        fn.kron_n_a(2 ** 2, np.transpose(spin2_i_z_t)))

    spin2_s_p_tl = np.kron(spin2_s_p_t, np.transpose(spin2_s_m_t)) - identity_mat + 0.5 * (
        fn.kron_a_n(spin2_s_z_t, 2 ** 2) +
        fn.kron_n_a(2 ** 2, np.transpose(spin2_s_z_t)))

    spin2_s_m_tl = np.kron(spin2_s_m_t, np.transpose(spin2_s_p_t)) - identity_mat - 0.5 * (
        fn.kron_a_n(spin2_s_z_t, 2 ** 2) +
        fn.kron_n_a(2 ** 2, np.transpose(spin2_s_z_t)))

    # Calculate relaxation matrices
    relax_t2_elec = (1 / param.t2_elec) * (np.kron(spin2_s_z_t, np.transpose(spin2_s_z_t)) - 0.5 * identity_mat)
    relax_t2_nuc = (1 / param.t2_nuc) * (np.kron(spin2_i_z_t, np.transpose(spin2_i_z_t)) - 0.5 * identity_mat)
    relax_t1 = gep * spin2_s_p_tl + gem * spin2_s_m_tl + gnp * spin2_i_p_tl + gnm * spin2_i_m_tl
    relax_mat = relax_t2_elec + relax_t2_nuc + relax_t1

    return relax_mat


def calculate_relaxation_mat_mance(eigvectors, eigvectors_inv, gnp, gnm, gep, gem, spin2_all):
    """ Calculate time dependent Liouville space relaxation matrix from Mance paper.
    """

    # Transform spin matrices into time dependent Hilbert space basis
    spin2_s_z_t, spin2_s_x_t,  spin2_s_p_t, spin2_s_m_t, spin2_i_z_t, spin2_i_x_t, spin2_i_p_t, spin2_i_m_t = \
        fn.basis_transform(eigvectors, eigvectors_inv, spin2_all)

    # Calculate Boltzmann factors and pre-allocate variables
    boltzmann_elec = param.electron_frequency * (sc.Planck / (sc.Boltzmann * param.temperature))
    boltzmann_nuc = param.freq_nuclear_1 * (sc.Planck / (sc.Boltzmann * param.temperature))
    relax_values_t1 = np.zeros((2 ** 2, 2 ** 2))
    relax_mat_t1 = np.zeros((4 ** 2, 4 ** 2))

    # Calculate T1 relaxation matrix using Boltzmann factors
    for a in range(0, 4):
        for b in range(0, 4):

            diff_elec = spin2_s_z_t[b, b] - spin2_s_z_t[a, a]
            diff_nuc = spin2_i_z_t[b, b] - spin2_i_z_t[a, a]

            boltzmann_factor = np.exp(-diff_elec * boltzmann_elec / 2 - diff_nuc * boltzmann_nuc / 2) / \
                               (np.exp(diff_elec * boltzmann_elec / 2 + diff_nuc * boltzmann_nuc / 2) +
                                np.exp(-diff_elec * boltzmann_elec / 2 - diff_nuc * boltzmann_nuc / 2))

            relax_values_t1[b, a] = (1 / param.t1_elec) * (spin2_s_x_t[b, a] * spin2_s_x_t[a, b] +
                                                           spin2_s_z_t[b, a] * spin2_s_z_t[a, b]) + \
                                    (1 / param.t1_nuc) * (spin2_i_x_t[b, a] * spin2_i_x_t[a, b] +
                                                          spin2_i_z_t[b, a] * spin2_i_z_t[a, b])

            relax_values_t1[b, a] = relax_values_t1[b, a] * boltzmann_factor

    # Set diagonals equal to zero
    for a in range(0, 4):
        relax_values_t1[a, a] = 0

    # Set diagonal elements
    for a in range(0, 4):
        for b in range(0, 4):
            if abs(a - b) > 0:
                relax_values_t1[b, b] = relax_values_t1[b, b] - relax_values_t1[a, b]

    # Transform to 16 by 16 matrix
    for a in range(0, 4):
        for b in range(0, 4):
            c = a * 4 + a
            d = b * 4 + b
            relax_mat_t1[d, c] = relax_values_t1[b, a]

    # Calculate T2 relaxation matrix using circular shift of spin matrices
    relax_values_t2 = np.zeros(4)
    for a in range(1, 4):
        relax_values_t2[a] = ((abs(spin2_s_z_t[a, a] - spin2_s_z_t[a - 1, a - 1])) ** 2) * (1 / param.t2_elec) + \
                        ((abs(spin2_i_z_t[a, a] - spin2_i_z_t[a - 1, a - 1])) ** 2) * (1 / param.t2_nuc)
        
    # Transform to 16 by 16 matrix
    relax_mat_t2 = np.diagflat(np.concatenate((relax_values_t2, np.roll(relax_values_t2, 1),
                                                np.roll(relax_values_t2, 2), np.roll(relax_values_t2, 3))))

    # Relaxation matrix as sum of t1 and t2 matrices
    relax_mat = -1*relax_mat_t2 + relax_mat_t1

    return relax_mat


if __name__ == "__main__":
    main()
    print('End of program.')
    plt.show()
