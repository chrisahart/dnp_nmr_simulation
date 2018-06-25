import parameters as param
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

# Pre-allocate arrays and calculate constants
hamiltonian = np.zeros((int(param.time_step_num), 4, 4), dtype=np.complex)
prop = np.zeros((int(param.time_step_num), 4, 4), dtype=np.complex)
nsp = 2  # number of spins

# G-anisotropy TODO move into function file
gx = param.electron_frequency * (2.00614 / 2) + 18.76e6
gy = param.electron_frequency * (2.00194 / 2) + 92.4e6
gz = param.electron_frequency * (2.00988 / 2) + 18.2e6

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

    hamiltonian[count, :] = (ganisohamil-param.electron_frequency)*param.spin2_s_z+\
                            param.freq_nuclear_1*param.spin2_i_z+hyperfine_zz + hyperfine_zx;

# Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
eigvals, eigvectors = np.linalg.eig(hamiltonian)
energies = np.real(eigvals)
eigvectors_inv = np.linalg.inv(eigvectors)

hamiltonian_ideal = 2 * np.pi * param.electron_frequency * param.spin2_s_z + \
                    2 * np.pi * -param.freq_nuclear_1 * param.spin2_i_z

# Calculate initial Zeeman basis density matrix from Boltzmann factors
boltzmann_factors = np.zeros(len(hamiltonian_ideal))
eigvals_ideal = np.linalg.eigvals(hamiltonian_ideal)
energies_ideal = np.real(np.copy(eigvals_ideal))
for count3 in range(0, hamiltonian.shape[1]):
    boltzmann_factors[count3] = np.exp((-param.hbar * energies_ideal[count3]) / (param.boltzmann * param.temperature))
density_mat = (1 / np.sum(boltzmann_factors)) * np.diagflat(boltzmann_factors)

# Calculate microwave Hamiltonian
microwave_hamiltonian_0 = 2 * np.pi * param.microwave_amplitude * param.spin2_s_x

# Propagation
for count in range(0, int(param.time_step_num)):

    # Express microwave Hamiltonian in initial Hamiltonian eigenbasis
    microwave_hamiltonian = np.matmul(eigvectors_inv[count], np.matmul(microwave_hamiltonian_0, eigvectors[count]))

    # Calculate total Hamiltonian
    energies_square = np.diagflat(energies[count])
    total_hamiltonian = energies_square + microwave_hamiltonian

    # Transform basis
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
    LD = np.kron(eigvectors[count], eigvectors[count])
    LDinv = np.kron(eigvectors_inv[count], eigvectors_inv[count])
    Ltot = Lhamilt + 1j * Rtot
    prop[count, :] = LD * la.expm(-1j * Ltot * param.time_step) * LDinv

# ?
prop_accu = np.eye(16)
for count in range(0, param.time_step_num):
    prop_accu = prop_accu * prop[count, :]

nrot = np.round(40 * param.freq_rotor)
iznew = np.zeros(nrot)
sznew = np.zeros(nrot)
rho0t = density_mat
for count in range(0, nrot):
    iznew[count] = np.trace(rho0t * param.spin2_i_z)
    sznew[count] = np.trace(rho0t * param.spin2_s_z)
    Lrho0t = np.reshape(rho0t, [16, 1])
    rho0t = prop_accu * Lrho0t
    rho0t = np.reshape(rho0t, [4, 4])

plt.plot(iznew)
plt.show()
