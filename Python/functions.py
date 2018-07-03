import numpy as np
import parameters as param
import spin_matrices as sp
from scipy import constants as sc
from scipy import linalg as la


def liouville_propagator(num_spins, energies, eigvectors, eigvectors_inv,
                         microwave_hamiltonian_init, calculate_relaxation_mat, spin_all):
    """ Calculate Liouville space propagator.
     """

    # Pre-allocate propagator
    propagator = np.zeros((int(param.time_step_num), 4 ** num_spins, 4 ** num_spins), dtype=np.complex)

    # Calculate variables for Liouville space relaxation
    p_e = np.tanh(0.5 * param.electron_frequency * (sc.Planck / (sc.Boltzmann * param.temperature)))
    p_n = np.tanh(0.5 * param.freq_nuclear_1 * (sc.Planck / (sc.Boltzmann * param.temperature)))
    gnp = 0.5 * (1 - p_n) * (1 / (1 * param.t1_nuc))
    gnm = 0.5 * (1 + p_n) * (1 / (1 * param.t1_nuc))
    gep = 0.5 * (1 - p_e) * (1 / (1 * param.t1_elec))
    gem = 0.5 * (1 + p_e) * (1 / (1 * param.t1_elec))

    for count in range(0, int(param.time_step_num)):

        # Transform microwave Hamiltonian into time dependent basis
        microwave_hamiltonian = np.matmul(eigvectors_inv[count], np.matmul(microwave_hamiltonian_init,
                                                                           eigvectors[count]))

        # Calculate total Hamiltonian
        total_hamiltonian = np.diag(energies[count]) + microwave_hamiltonian

        # Transform Hilbert space Hamiltonian into Liouville space
        hamiltonian_liouville = kron_a_n(total_hamiltonian, 2 ** num_spins) - \
                                kron_n_a(2 ** num_spins, np.transpose(total_hamiltonian))

        # Calculate time dependent Liouville space relaxation matrix
        relax_mat = calculate_relaxation_mat(eigvectors[count], eigvectors_inv[count], gnp, gnm, gep, gem, spin_all)

        # Calculate Liouville space eigenvectors
        eigvectors_liouville = np.kron(eigvectors[count], eigvectors[count])
        eigvectors_inv_liouville = np.kron(eigvectors_inv[count], eigvectors_inv[count])

        # Calculate Liouville space propagator
        liouvillian = hamiltonian_liouville + 1j * relax_mat

        test2 = la.cosm(liouvillian * param.time_step) - 1j * la.sinm(liouvillian * param.time_step)
        propagator[count, :] = np.matmul(eigvectors_inv_liouville, np.matmul(test2, eigvectors_liouville))

        # propagator[count, :] = np.matmul(eigvectors_inv_liouville,
        #                                  np.matmul(la.expm(-1j * liouvillian * param.time_step), eigvectors_liouville))

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


def kron_a_n(A, N):  # Simulates np.kron(A, np.eye(N))

    m,n = A.shape
    out = np.zeros((m,N,n,N),dtype=A.dtype)
    r = np.arange(N)
    out[:,r,:,r] = A
    out.shape = (m*N,n*N)
    return out


def kron_n_a(N, A):  # Simulates np.kron(np.eye(N), A)
    m,n = A.shape
    out = np.zeros((N,m,N,n),dtype=A.dtype)
    r = np.arange(N)
    out[r,:,r,:] = A
    out.shape = (m*N,n*N)
    return out
