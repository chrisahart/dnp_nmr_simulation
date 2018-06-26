import parameters as param
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

# Pre-allocate arrays and calculate constants
hamiltonian = np.zeros((int(param.time_step_num), 4, 4), dtype=np.complex)
prop = np.zeros((int(param.time_step_num), 4**2, 4**2), dtype=np.complex)
nsp = 2  # number of spins

# G-anisotropy TODO move into function file
gx = param.electron_frequency * (2.00614 / 2) + 18.76e6
gy = param.electron_frequency * (2.00194 / 2) + 92.4e6
gz = param.electron_frequency * (2.00988 / 2) + 18.2e6

print(gx)

ca = np.cos(param.orientation_tempol[0])
cb = np.cos(param.orientation_tempol[1])
cg = np.cos(param.orientation_tempol[2])
sa = np.sin(param.orientation_tempol[0])
sb = np.sin(param.orientation_tempol[1])
sg = np.sin(param.orientation_tempol[2])

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

print(c0)
print(c1)
print(c2)
print(c3)
print(c4)

# Construct intrinsic Hamiltonian (Hilbert space)
for count in range(0, int(param.time_step_num)):

    # TODO Add all of below into function to create Hamiltonian
    hyperfine_zz = param.hyperfine_coupling * (-0.5 * (np.sin(param.hyperfine_angles_1[1]) ** 2) * np.cos(
        2 * (2 * np.pi * param.freq_rotor * count * param.time_step + param.hyperfine_angles_1[2])) + np.sqrt(
        2) * np.sin(param.hyperfine_angles_1[1]) * np.cos(param.hyperfine_angles_1[1]) * np.cos(
        2 * np.pi * param.freq_rotor * count * param.time_step + param.hyperfine_angles_1[2])) * 2 * param.IzSz

    hyperfine_zx = param.hyperfine_coupling * (
        -0.5 * np.sin(param.hyperfine_angles_1[1]) * np.cos(param.hyperfine_angles_1[1]) * np.cos(
            2 * np.pi * param.freq_rotor * count * param.time_step + param.hyperfine_angles_1[2]) - (np.sqrt(2) / 4) * (
            np.sin(param.hyperfine_angles_1[1]) ** 2) * np.cos(
            2 * (2 * np.pi * param.freq_rotor * count * param.time_step + param.hyperfine_angles_1[2])) + (
            np.sqrt(2) / 4) * (3 * (np.cos(param.hyperfine_angles_1[1]) ** 2) - 1)) * (param.IpSz + param.ImSz)

    ganisohamil = c0 + c1 * np.cos(2 * np.pi * param.freq_rotor * count * param.time_step) + c2 * np.sin(
        2 * np.pi * param.freq_rotor * count * param.time_step) + c3 * np.cos(
        2 * np.pi * param.freq_rotor * count * param.time_step * 2) + c4 * np.sin(
        2 * np.pi * param.freq_rotor * count * param.time_step * 2)

    # hamiltonian[count, :] = (ganisohamil-param.wme)*param.spin2_s_z

    hamiltonian[count, :] = (ganisohamil - param.wme) * param.spin2_s_z + \
                            param.freq_nuclear_1 * param.spin2_i_z #+ hyperfine_zz + hyperfine_zx

# Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
eigvals, eigvectors = np.linalg.eig(hamiltonian)
energies = np.real(eigvals)
eigvectors_inv = np.linalg.inv(eigvectors)

print('(ganisohamil) \n', (ganisohamil))
print('(param.wme) \n', (param.wme))
print('(ganisohamil-param.electron_frequency) \n', (ganisohamil-param.wme))

print('hamiltonian \n', hamiltonian.shape)
print('hamiltonian[0,0] \n', hamiltonian[0, :, :])

print('energies \n', energies.shape)
print('energies[0,0] \n', energies[0,:])

print('energies \n', eigvectors.shape)
print('eigvectors[0,0] \n', eigvectors[0, :, :])

# hamiltonian_ideal = param.electron_frequency * param.spin2_s_z + \
#                     param.freq_nuclear_1 * param.spin2_i_z

# Calculate initial Zeeman basis density matrix from Boltzmann factors
# boltzmann_factors = np.zeros(len(hamiltonian_ideal))
# eigvals_ideal = np.linalg.eigvals(hamiltonian_ideal)
# energies_ideal = np.real(np.copy(eigvals_ideal))
# for count3 in range(0, hamiltonian.shape[1]):
#     boltzmann_factors[count3] = np.exp((-param.hbar * 2 * np.pi * energies_ideal[count3]) / (param.boltzmann * param.temperature))
# density_mat = (1 / np.sum(boltzmann_factors)) * np.diagflat(boltzmann_factors)
# print('np.sum(boltzmann_factors) \n', np.sum(boltzmann_factors))

planck = 2 * np.pi * param.hbar

Zp = np.exp(
    (param.electron_frequency + param.freq_nuclear_1) * planck / (param.boltzmann * param.temperature)) + np.exp(
    (-param.electron_frequency + param.freq_nuclear_1) * planck / (param.boltzmann * param.temperature)) + np.exp(
    (param.electron_frequency - param.freq_nuclear_1) * planck / (param.boltzmann * param.temperature)) + np.exp(
    -(param.electron_frequency + param.freq_nuclear_1) * planck / (param.boltzmann * param.temperature))
rho_zeeman = (1 / Zp) * la.expm(
    -(param.electron_frequency * param.spin2_s_z + param.freq_nuclear_1 * param.spin2_i_z) * planck / (
        param.boltzmann * param.temperature))
density_mat = rho_zeeman

print('density_mat \n', density_mat)

# Calculate microwave Hamiltonian
microwave_hamiltonian_0 = param.microwave_amplitude * param.spin2_s_x
print('microwave_hamiltonian_0 \n', microwave_hamiltonian_0)

# Propagation
#for count in range(0, int(param.time_step_num)):
for count in range(0, 1):

    # Express microwave Hamiltonian in initial Hamiltonian eigenbasis
    microwave_hamiltonian = np.matmul(eigvectors_inv[count], np.matmul(microwave_hamiltonian_0, eigvectors[count]))
    print('microwave_hamiltonian \n', microwave_hamiltonian)

    # Calculate total Hamiltonian
    energies_square = np.diagflat(energies[count])
    total_hamiltonian = energies_square + microwave_hamiltonian

    # Transform basis TODO move into function file def(transform_spin_matrices)
    Sxt = np.matmul(eigvectors_inv[count], np.matmul(param.spin2_s_x, eigvectors[count]))
    Syt = np.matmul(eigvectors_inv[count], np.matmul(param.spin2_s_y, eigvectors[count]))
    Szt = np.matmul(eigvectors_inv[count], np.matmul(param.spin2_s_z, eigvectors[count]))

    # Transform basis
    Ixt = np.matmul(eigvectors_inv[count], np.matmul(param.spin2_i_x, eigvectors[count]))
    Iyt = np.matmul(eigvectors_inv[count], np.matmul(param.spin2_i_y, eigvectors[count]))
    Izt = np.matmul(eigvectors_inv[count], np.matmul(param.spin2_i_z, eigvectors[count]))

    # Transform basis
    Spt = np.matmul(eigvectors_inv[count], np.matmul(param.Sp, eigvectors[count]))
    Smt = np.matmul(eigvectors_inv[count], np.matmul(param.Sm, eigvectors[count]))
    Ipt = np.matmul(eigvectors_inv[count], np.matmul(param.Ip, eigvectors[count]))
    Imt = np.matmul(eigvectors_inv[count], np.matmul(param.Im, eigvectors[count]))

    # Liouvillian?
    LSzt = np.kron(Szt, np.transpose(Szt))
    LIzt = np.kron(Izt, np.transpose(Izt))
    RSzt = np.kron(np.eye(2 ** nsp), np.eye(2 ** nsp))
    RIzt = np.kron(np.eye(2 ** nsp), np.eye(2 ** nsp))

    # Liouvillian?
    LIpt = 1 * np.kron(Ipt, np.transpose(Imt)) - 0.5 * np.eye(4 ** nsp) + .5 * (
        np.kron(Izt, np.eye(2 ** nsp)) + np.kron(np.eye(2 ** nsp), np.transpose(Izt)))
    LImt = 1 * np.kron(Imt, np.transpose(Ipt)) - 0.5 * np.eye(4 ** nsp) + .5 * (
        np.kron(Izt, np.eye(2 ** nsp)) + np.kron(np.eye(2 ** nsp), np.transpose(Izt)))
    LSpt = 1 * np.kron(Spt, np.transpose(Smt)) - 0.5 * np.eye(4 ** nsp) + .5 * (
        np.kron(Szt, np.eye(2 ** nsp)) + np.kron(np.eye(2 ** nsp), np.transpose(Szt)))
    LSmt = 1 * np.kron(Smt, np.transpose(Spt)) - 0.5 * np.eye(4 ** nsp) + .5 * (
        np.kron(Szt, np.eye(2 ** nsp)) + np.kron(np.eye(2 ** nsp), np.transpose(Szt)))

    # Liouvillian?
    Rt2ematv3 = (1 / param.t2_elec) * (LSzt - 0.25 * RSzt)
    Rt2nmatv3 = (1 / param.t2_nuc) * (LIzt - 0.25 * RIzt)
    Rt1mat = param.gep * LSpt + param.gem * LSmt + param.gnp * LIpt + param.gnm * LImt
    Rtot = 1 * Rt1mat + 1 * Rt2ematv3 + 1 * Rt2nmatv3

    # Liouvillian?
    Lhamilt = np.kron(total_hamiltonian, np.eye(4)) - np.kron(np.eye(4), np.transpose(total_hamiltonian))
    print('Lhamilt \n', Lhamilt.shape)
    print('Lhamilt \n', Lhamilt)

    LD = np.kron(eigvectors[count], eigvectors[count])
    LDinv = np.kron(eigvectors_inv[count], eigvectors_inv[count])
    Ltot = Lhamilt + 1j * Rtot
    test = LD * la.expm(-1j * Ltot * param.time_step) * LDinv
    prop[count, :] = LD * la.expm(-1j * Ltot * param.time_step) * LDinv
    # print('test \n', test.shape)
    # print('test \n', test)

# # ?
# prop_accu = np.eye(16)
# for count in range(0, int(param.time_step_num)):
#     prop_accu = prop_accu * prop[count, :]
#     # print('prop[count, :] \n', prop[count, :])
#     # print('prop_accu \n', prop_accu)
#
# nrot = int(np.round(40 * param.freq_rotor))
# iznew = np.zeros(nrot)
# sznew = np.zeros(nrot)
# rho0t = density_mat
#
# for count in range(0, 1):
#     iznew[count] = np.trace(rho0t * param.spin2_i_z)
#     sznew[count] = np.trace(rho0t * param.spin2_s_z)
#     Lrho0t = np.reshape(rho0t, [16, 1])
#
#     print('iznew \n', iznew[count])  # agrees with Matlab
#     print('sznew \n', sznew[count])  # agrees with Matlab
#
#     print('Lrho0t \n', Lrho0t.shape)  # agrees with Matlab
#     print('Lrho0t \n', Lrho0t)  # agrees with Matlab
#
#     rho0t = np.matmul(prop_accu, Lrho0t)
#
#     # print('prop_accu \n', prop_accu.shape)  # does not agree with matlab
#     # print('prop_accu[0, 0] \n', prop_accu[0, 0])
#     # print('prop_accu \n', prop_accu[1, 1])
#     # print('prop_accu \n', prop_accu)
#
#     # print('rho0t after \n', rho0t.shape)
#     # print('rho0t after \n', rho0t)
#
#     rho0t = np.reshape(rho0t, [4, 4])
#     print('rho0t final \n', rho0t.shape)
#     print('rho0t final \n', rho0t)
#
# plt.plot(iznew)
# plt.show()
