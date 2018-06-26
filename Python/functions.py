import numpy as np
import parameters as param
import spin_matrices as sp
from scipy import constants as sc
from scipy import linalg as la


def density_mat_thermal():
    """ Calculate initial thermal density matrix from Boltzmann factors.
    """

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

    return density_mat


def anisotropy_coefficients(angles):
    """ Calculate time independent coefficients for electron g-anisotropy.
    """

    gx = param.electron_frequency * (2.00614 / 2) + 18.76e6
    gy = param.electron_frequency * (2.00194 / 2) + 92.4e6
    gz = param.electron_frequency * (2.00988 / 2) + 18.2e6

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

    ganisohamil = c0 + c1 * np.cos(2 * np.pi * param.freq_rotor * time) + c2 * np.sin(
        2 * np.pi * param.freq_rotor * time) + c3 * np.cos(
        2 * np.pi * param.freq_rotor * time * 2) + c4 * np.sin(
        2 * np.pi * param.freq_rotor * time * 2)

    return ganisohamil


def hyperfine(time):
    """ Calculate time dependent hyperfine, neglecting zy interaction.
    """

    hyperfine_zz = param.hyperfine_coupling * (-0.5 * (np.sin(param.hyperfine_angles_1[1]) ** 2) * np.cos(
        2 * (2 * np.pi * param.freq_rotor * time + param.hyperfine_angles_1[2])) + np.sqrt(
        2) * np.sin(param.hyperfine_angles_1[1]) * np.cos(param.hyperfine_angles_1[1]) * np.cos(
        2 * np.pi * param.freq_rotor * time + param.hyperfine_angles_1[2])) * 2 * np.matmul(sp.spin2_i_z,
                                                                                            sp.spin2_s_z)

    hyperfine_zx = param.hyperfine_coupling * (
            -0.5 * np.sin(param.hyperfine_angles_1[1]) * np.cos(param.hyperfine_angles_1[1]) * np.cos(
        2 * np.pi * param.freq_rotor * time + param.hyperfine_angles_1[2]) - (np.sqrt(2) / 4) * (
                    np.sin(param.hyperfine_angles_1[1]) ** 2) * np.cos(
        2 * (2 * np.pi * param.freq_rotor * time + param.hyperfine_angles_1[2])) + (
                    np.sqrt(2) / 4) * (3 * (np.cos(param.hyperfine_angles_1[1]) ** 2) - 1)) * \
                   (np.matmul(sp.spin2_i_p, sp.spin2_s_z) + np.matmul(sp.spin2_i_m, sp.spin2_s_z))

    return hyperfine_zz, hyperfine_zx
