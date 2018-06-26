import parameters as param
import numpy as np
from scipy import linalg as la
from scipy import constants as sc
import spin_matrices as sp
import functions as fn


def dynamics():
    """ Construct Hamiltonian.
    """

    # Calculate variables
    time_end = 1 / param.freq_rotor
    time_step = time_end / param.time_step_num
    c0, c1, c2, c3, c4 = fn.anisotropy_coefficients(param.orientation_tempol)

    # Pre-allocate arrays
    hamiltonian = np.zeros((int(param.time_step_num), 4, 4), dtype=np.complex)
    prop = np.zeros((int(param.time_step_num), 4**2, 4**2), dtype=np.complex)

    # Construct intrinsic Hamiltonian (Hilbert space)
    for count in range(0, int(param.time_step_num)):

        # TODO Add all of below into function to create Hamiltonian
        hyperfine_zz = param.hyperfine_coupling * (-0.5 * (np.sin(param.hyperfine_angles_1[1]) ** 2) * np.cos(
            2 * (2 * np.pi * param.freq_rotor * count * time_step + param.hyperfine_angles_1[2])) + np.sqrt(
            2) * np.sin(param.hyperfine_angles_1[1]) * np.cos(param.hyperfine_angles_1[1]) * np.cos(
            2 * np.pi * param.freq_rotor * count * time_step + param.hyperfine_angles_1[2])) * 2 * np.matmul(sp.spin2_i_z, sp.spin2_s_z)

        hyperfine_zx = param.hyperfine_coupling * (
            -0.5 * np.sin(param.hyperfine_angles_1[1]) * np.cos(param.hyperfine_angles_1[1]) * np.cos(
                2 * np.pi * param.freq_rotor * count * time_step + param.hyperfine_angles_1[2]) - (np.sqrt(2) / 4) * (
                np.sin(param.hyperfine_angles_1[1]) ** 2) * np.cos(
                2 * (2 * np.pi * param.freq_rotor * count * time_step + param.hyperfine_angles_1[2])) + (
                np.sqrt(2) / 4) * (3 * (np.cos(param.hyperfine_angles_1[1]) ** 2) - 1)) * \
                       (np.matmul(sp.spin2_i_p, sp.spin2_s_z) + np.matmul(sp.spin2_i_m, sp.spin2_s_z))

        ganisohamil = c0 + c1 * np.cos(2 * np.pi * param.freq_rotor * count * time_step) + c2 * np.sin(
            2 * np.pi * param.freq_rotor * count * time_step) + c3 * np.cos(
            2 * np.pi * param.freq_rotor * count * time_step * 2) + c4 * np.sin(
            2 * np.pi * param.freq_rotor * count * time_step * 2)

        hamiltonian[count, :] = (ganisohamil - param.microwave_frequency) * sp.spin2_s_z + \
                                param.freq_nuclear_1 * sp.spin2_i_z + hyperfine_zz + hyperfine_zx

    # Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
    eigvals, eigvectors = np.linalg.eig(hamiltonian)
    energies = np.real(eigvals)
    eigvectors_inv = np.linalg.inv(eigvectors)

    # hamiltonian_ideal = param.electron_frequency * sp.spin2_i_z + \
    #                     param.freq_nuclear_1 * sp.spin2_i_z

    # Calculate initial Zeeman basis density matrix from Boltzmann factors
    # boltzmann_factors = np.zeros(len(hamiltonian_ideal))
    # eigvals_ideal = np.linalg.eigvals(hamiltonian_ideal)
    # energies_ideal = np.real(np.copy(eigvals_ideal))
    # for count3 in range(0, hamiltonian.shape[1]):
    #     boltzmann_factors[count3] = np.exp((-param.sc.hbar * 2 * np.pi * energies_ideal[count3]) / (sc.Boltzmann * param.temperature))
    # density_mat = (1 / np.sum(boltzmann_factors)) * np.diagflat(boltzmann_factors)
    # print('np.sum(boltzmann_factors) \n', np.sum(boltzmann_factors))

    Zp = np.exp(
        (param.electron_frequency + param.freq_nuclear_1) * sc.Planck / (sc.Boltzmann * param.temperature)) + np.exp(
        (-param.electron_frequency + param.freq_nuclear_1) * sc.Planck / (sc.Boltzmann * param.temperature)) + np.exp(
        (param.electron_frequency - param.freq_nuclear_1) * sc.Planck / (sc.Boltzmann * param.temperature)) + np.exp(
        -(param.electron_frequency + param.freq_nuclear_1) * sc.Planck / (sc.Boltzmann * param.temperature))
    rho_zeeman = (1 / Zp) * la.expm(
        -(param.electron_frequency * sp.spin2_s_z + param.freq_nuclear_1 * sp.spin2_i_z) * sc.Planck / (
            sc.Boltzmann * param.temperature))
    density_mat = rho_zeeman

    # Calculate microwave Hamiltonian
    microwave_hamiltonian_0 = param.microwave_amplitude * sp.spin2_s_z

    # Propagation
    for count in range(0, int(param.time_step_num)):

        # Express microwave Hamiltonian in initial Hamiltonian eigenbasis
        microwave_hamiltonian = np.matmul(eigvectors_inv[count], np.matmul(microwave_hamiltonian_0, eigvectors[count]))

        # Calculate total Hamiltonian
        energies_square = np.diagflat(energies[count])
        total_hamiltonian = energies_square + microwave_hamiltonian

        # Transform basis TODO move into function file def(transform_spin_matrices)
        Sxt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_s_x, eigvectors[count]))
        Syt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_s_y, eigvectors[count]))
        Szt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_s_z, eigvectors[count]))

        #print('Szt \n', Szt)

        # Transform basis
        Ixt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_i_x, eigvectors[count]))
        Iyt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_i_y, eigvectors[count]))
        Izt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_i_z, eigvectors[count]))

        # Transform basis
        Spt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_s_p, eigvectors[count]))
        Smt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_s_m, eigvectors[count]))
        Ipt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_i_p, eigvectors[count]))
        Imt = np.matmul(eigvectors_inv[count], np.matmul(sp.spin2_i_m, eigvectors[count]))

        # Liouvillian?
        LSzt = np.kron(Szt, np.transpose(Szt))
        LIzt = np.kron(Izt, np.transpose(Izt))
        RSzt = np.kron(np.eye(2 ** param.num_spins), np.eye(2 ** param.num_spins))
        RIzt = np.kron(np.eye(2 ** param.num_spins), np.eye(2 ** param.num_spins))

        # Liouvillian?
        LIpt = np.kron(Ipt, np.transpose(Imt)) - 0.5 * np.eye(4 ** param.num_spins) + 0.5 * (
            np.kron(Izt, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins), np.transpose(Izt)))

        LImt = np.kron(Imt, np.transpose(Ipt)) - 0.5 * np.eye(4 ** param.num_spins) - 0.5 * (
            np.kron(Izt, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins), np.transpose(Izt)))
        LSpt = np.kron(Spt, np.transpose(Smt)) - 0.5 * np.eye(4 ** param.num_spins) + 0.5 * (
            np.kron(Szt, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins), np.transpose(Szt)))
        LSmt = np.kron(Smt, np.transpose(Spt)) - 0.5 * np.eye(4 ** param.num_spins) - 0.5 * (
            np.kron(Szt, np.eye(2 ** param.num_spins)) + np.kron(np.eye(2 ** param.num_spins), np.transpose(Szt)))

        # Liouvillian?
        Rt2ematv3 = (1 / param.t2_elec) * (LSzt - 0.25 * RSzt)
        Rt2nmatv3 = (1 / param.t2_nuc) * (LIzt - 0.25 * RIzt)
        Rt1mat = param.gep * LSpt + param.gem * LSmt + param.gnp * LIpt + param.gnm * LImt
        Rtot = Rt1mat + Rt2ematv3 + Rt2nmatv3

        # Liouvillian?
        Lhamilt = np.kron(total_hamiltonian, np.eye(4)) - np.kron(np.eye(4), np.transpose(total_hamiltonian))

        LD = np.kron(eigvectors[count], eigvectors[count])
        LDinv = np.kron(eigvectors_inv[count], eigvectors_inv[count])
        Ltot = Lhamilt + 1j * Rtot
        prop[count, :] = np.matmul(LDinv, np.matmul(la.expm(-1j * Ltot * time_step), LD))

    prop_accu = np.eye(16)
    for count in range(0, int(param.time_step_num)):
        prop_accu = np.matmul(prop_accu, prop[count, :])

    nrot = int(np.round(40 * param.freq_rotor))
    iznew = np.zeros(nrot)
    sznew = np.zeros(nrot)
    rho0t = density_mat

    for count in range(0, nrot):
        iznew[count] = np.trace(np.real(rho0t) * sp.spin2_i_z)
        sznew[count] = np.trace(np.real(rho0t) * sp.spin2_s_z)
        Lrho0t = np.reshape(rho0t, [16, 1])

        rho0t = np.matmul(prop_accu, Lrho0t)
        rho0t = np.reshape(rho0t, [4, 4])

    return energies, iznew, sznew
