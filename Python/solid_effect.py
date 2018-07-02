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
        Can't use Numba or Pythran as not pure python, using Numpy instead.
        Multiprocessing doesn't work due to overhead, and nature of numpy matrix products.
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
    """ Calculate dynamics of e-n system, calling required functions.
    """

    # Pre-allocate arrays
    propagator_strobe = np.eye(4 ** 2)

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
    # hamiltonian = fortran.f2py_dynamics.spin2_electronuclear()

    # Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
    eigvals, eigvectors = np.linalg.eig(hamiltonian)
    energies = np.real(eigvals)
    eigvectors_inv = np.linalg.inv(eigvectors)

    # Calculate microwave Hamiltonian
    microwave_hamiltonian_init = microwave_amplitude * spin2_s_x

    # Calculate Liouville space propagator with relaxation
    propagator = fn.liouville_propagator(2, energies, eigvectors, eigvectors_inv,
                                         microwave_hamiltonian_init, calculate_relaxation_mat, spin2_all)

    # Propagate density matrix for single rotor period, calculating polarisations
    pol_i_z_rot, pol_s_z_rot = calculate_polarisation_sub_rotor(density_mat, propagator, spin2_s_z, spin2_i_z)

    # Calculate stroboscopic propagator (product of all operators within rotor period)
    for count in range(0, int(param.time_step_num)):
        propagator_strobe = np.matmul(propagator_strobe, propagator[count, :])

    # Propagate density matrix stroboscopically, calculating polarisations
    pol_i_z, pol_s_z = calculate_polarisation_rotor(density_mat, propagator_strobe, spin2_s_z, spin2_i_z)

    return pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies


def calculate_hamiltonian(spin2_s_z, spin2_i_x, spin2_i_z):
    """ Calculate Hamiltonian of e-n system.
    """

    # Pre-allocate hamiltonian
    hamiltonian = np.zeros((int(param.time_step_num), 2 ** 2, 2 ** 2,), dtype=np.complex)

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
    spin2_s_z_t, spin2_s_p_t, spin2_s_m_t, spin2_i_z_t, spin2_i_p_t, spin2_i_m_t = \
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


if __name__ == "__main__":
    main()
    print('End of program.')
    plt.show()
