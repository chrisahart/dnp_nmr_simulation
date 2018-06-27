import parameters as param
import numpy as np
from scipy import linalg as la
from scipy import constants as sc
import spin_matrices as sp
import functions as fn


def dynamics(microwave_amplitude):
    """ Solid effect dynamics, needs generalizing and optimising.
    """

    # Pre-allocate arrays for Hamiltonian and propagator
    hamiltonian = np.zeros((int(param.time_step_num), param.num_spins ** 2, param.num_spins ** 2), dtype=np.complex)
    prop = np.zeros((int(param.time_step_num), param.num_spins ** 4, param.num_spins ** 4), dtype=np.complex)

    # Calculate time independent electron g-anisotropy coefficients
    c0, c1, c2, c3, c4 = fn.anisotropy_coefficients(param.orientation_tempol)  # TODO generalize

    # Calculate thermal density matrix
    density_mat = fn.density_mat_thermal()  # TODO generalize

    # Construct intrinsic Hamiltonian (Hilbert space)
    for count in range(0, int(param.time_step_num)):

        hyperfine_zz, hyperfine_zx = fn.hyperfine(count * param.time_step)  # TODO generalize

        ganisotropy = fn.anisotropy(c0, c1, c2, c3, c4, count * param.time_step)  # TODO generalize

        hamiltonian[count, :] = (ganisotropy - param.microwave_frequency) * sp.spin2_s_z + \
                                param.freq_nuclear_1 * sp.spin2_i_z + hyperfine_zz + hyperfine_zx  # TODO generalize

    # Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
    eigvals, eigvectors = np.linalg.eig(hamiltonian)
    energies = np.real(eigvals)
    eigvectors_inv = np.linalg.inv(eigvectors)

    # TODO generalize (add to function to sum to number of electrons)
    microwave_hamiltonian_0 = microwave_amplitude * sp.spin2_s_x

    # Propagate density matrix
    for count in range(0, int(param.time_step_num)):

        microwave_hamiltonian = np.matmul(eigvectors_inv[count],
                                          np.matmul(microwave_hamiltonian_0, eigvectors[count]))

        # Calculate total Hamiltonian
        total_hamiltonian = np.diagflat(energies[count]) + microwave_hamiltonian

        # Transform spin matrices into time dependent basis TODO generalize
        Sxt, Syt, Szt, Ixt, Iyt, Izt, Spt, Smt, Ipt, Imt = fn.basis_transform(eigvectors, eigvectors_inv, count)

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
        prop[count, :] = np.matmul(LDinv, np.matmul(la.expm(-1j * Ltot * param.time_step), LD))

    prop_accu = np.eye(16)
    for count in range(0, int(param.time_step_num)):
        prop_accu = np.matmul(prop_accu, prop[count, :])

    iznew = np.zeros(param.nrot)
    sznew = np.zeros(param.nrot)
    rho0t = density_mat

    for count in range(0, param.nrot):
        iznew[count] = np.trace(np.real(rho0t) * sp.spin2_i_z)
        sznew[count] = np.trace(np.real(rho0t) * sp.spin2_s_z)
        Lrho0t = np.reshape(rho0t, [16, 1])

        rho0t = np.matmul(prop_accu, Lrho0t)
        rho0t = np.reshape(rho0t, [4, 4])

    return iznew, sznew
