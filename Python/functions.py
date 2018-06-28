import numpy as np
import parameters as param
import spin_matrices as sp
from scipy import constants as sc
from scipy import linalg as la


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

    gx = param.electron_frequency * param.gtensor[0] + param.nitrogen_coupling[0]
    gy = param.electron_frequency * param.gtensor[1] + param.nitrogen_coupling[1]
    gz = param.electron_frequency * param.gtensor[2] + param.nitrogen_coupling[2]

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
                                               np.cos(2 * np.pi * param.freq_rotor * time + hyperfine_angles[2])) * 2

    hyperfine_zx = param.hyperfine_coupling * (-0.5 * np.sin(hyperfine_angles[1]) * np.cos(hyperfine_angles[1]) *
                                               np.cos(2 * np.pi * param.freq_rotor * time + hyperfine_angles[2]) -
                                               (np.sqrt(2) / 4) * (np.sin(hyperfine_angles[1]) ** 2) *
                                               np.cos(2 * (2 * np.pi * param.freq_rotor * time + hyperfine_angles[2])) +
                                               (np.sqrt(2) / 4) * (3 * (np.cos(hyperfine_angles[1]) ** 2) - 1))

    return hyperfine_zz, hyperfine_zx
