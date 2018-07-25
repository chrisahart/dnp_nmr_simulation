import numpy as np
from scipy import constants as sc
import spin_matrices as sp

# Literature parameters
b_field = 9.4                                                       # Magnetic field
temperature = 100                                                   # Temperature
freq_nuclear_1 = -400.9E6                                           # todo change to nuclear_frequency
electron_frequency = 28.025E9 * b_field                             # Electron frequency
microwave_frequency = 264E9 #265.2E9 #                                         # Microwave frequency
freq_rotor = 3E3                                                    # Rotor frequency
orientation_se = np.radians((253.6, 105.1, 123.8))                  # G anisotropy angles for electron 1 (SE)
orientation_ce_1 = np.radians((253.6, 105.1, 123.8))                # G anisotropy angles for electron 1 (CE)
orientation_ce_2 = orientation_ce_1 + np.radians((102, 104, 124))   # G anisotropy angles for electron 2 (CE)
gtensor = np.array([(2.00614 / 2), (2.00194 / 2), (2.00988 / 2)])   # G tensor principal values
hyperfine_coupling = 3E6                                            # Hyperfine coupling amplitude
hyperfine_angles_1 = np.radians((0, 0, 0))                          # Hyperfine angles for e1-n
#hyperfine_angles_1 = np.radians((253.6, 105.1, 123.8))              # Hyperfine angles for e1-n
t1_elec = 0.3E-3                                                    # Electron spin-lattice relaxation T1 (s)
t1_nuc = 10                                                         # Nuclear spin-lattice relaxation T1 (s)
t2_elec = 1e-6                                                      # T2 electron
t2_nuc = 1e-3                                                       # T2 nucleus
microwave_amplitude = np.ones(2) * 0.85 * 1E6 #np.array([0.85]) * 1E6                        # Microwave field amplitude np.arange(1, 10, 5)
#microwave_amplitude = np.arange(1, 200, 5) * 1E6
#microwave_amplitude = np.array([75]) * 1E6

# System variables
time_step_num = int(1E4)                                            # Number of timesteps within rotor period
time_step = (1 / freq_rotor) / time_step_num                        # Value of timestep within rotor period
num_periods = 40                                                    # Time in seconds to propagate system  todo rename
num_timesteps_prop = int(np.round(num_periods * freq_rotor))        # Number of timesteps to propagate system
