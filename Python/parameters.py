import numpy as np

# Universal constants
boltzmann = 1.38E-23                                                # Boltzmann constant
avagadro = 6.022E23                                                 # Avogadro's constant
hbar = 1.054E-34                                                    # Reduced Planck constant

# Literature parameters
b_field = 9.4                                                       # Magnetic field
temperature = 25                                                           # Temperature
freq_nuclear_1 = -400.9E6                                           # Nuclear frequency 1
electron_frequency = 28.025E9 * b_field                             # Electron frequency
microwave_frequency = 264E9                                         # Microwave frequency
freq_rotor = 3E3                                                    # Rotor frequency
microwave_amplitude = 1E6                                        # Microwave field amplitude
orientation_tempol = np.radians((253.6, 105.1, 123.8))                    # G anisotropy angles for Tempol (SE)
gtensor = np.array([(2.00614 / 2), (2.00194 / 2), (2.00988 / 2)])  # G tensor principal values
hyperfine_coupling = 3E6                                            # Hyperfine coupling amplitude
hyperfine_angles_1 = np.radians([0, 0, 0])                          # Hyperfine angles for e1-n
t1_elec = 0.3E-3;                       # Electron spin-lattice relaxtion T1 (s)
t1_nuc = 10;                            # Nuclear spin-lattice relaxation T1 (s)
t2_elec=1e-6; # T2 electron
t2_nuc=1e-3; # T2 nucleus

# What are these, something to do with Liouvillian?
t_corr = hbar / (boltzmann * temperature)
p_e = np.tanh(0.5 * electron_frequency * t_corr)
p_n = np.tanh(0.5 * freq_nuclear_1 * t_corr)
gnp = 0.5 * (1 - p_n) * (1 / (1 * t1_nuc))
gnm = 0.5 * (1 + p_n) * (1 / (1 * t1_nuc))
gep = 0.5 * (1 - p_e) * (1 / (1 * t1_elec))
gem = 0.5 * (1 + p_e) * (1 / (1 * t1_elec))

# System variables
time_step_num = 1E4
time_end = 1/freq_rotor
time_step = time_end/time_step_num

# 2x2 Matrix operators for spin ensemble TODO move to matrices file, is there a way to automate construction?
spin1_x = 1/2 * np.array([[0, 1],
                          [1, 0]])

spin1_y = 1/2 * np.array([[0, -1j],
                          [1j, 0]])

spin1_z = 1/2 * np.array([[1, 0],
                          [0, -1]])

# 4x4 Matrix operators for S operator
spin2_s_x = np.kron(spin1_x, np.identity(2))
spin2_s_y = np.kron(spin1_y, np.identity(2))
spin2_s_z = np.kron(spin1_z, np.identity(2))

# 4x4 Matrix operators for I operator
spin2_i_x = np.kron(np.identity(2), spin1_x)
spin2_i_y = np.kron(np.identity(2), spin1_y)
spin2_i_z = np.kron(np.identity(2), spin1_z)

# TODO Add + and - spin operators
Sp = 1
Sm = 1
Ip = 1
Im = 1
# Sp=sop(spins,'+e');
# Sm=sop(spins,'-e');
# Ip=sop(spins,'e+');
# Im=sop(spins,'e-');

# TODO Add
IzSz = np.matmul(spin2_i_z, spin2_s_z)
IpSz = np.matmul(Ip, spin2_s_z)
ImSz = np.matmul(Im, spin2_s_z)
