import numpy as np
import matplotlib.pyplot as plt
import solid_effect
import time
import parameters as param
import solid_effect_plotting
from shutil import copyfile
import os

""" Call dynamics function as required, saving output to disk.
    Effect of microwave amplitude on nuclear enhancement can be simulated.
"""

# Pre-allocate arrays
pol_nuc = np.zeros((param.microwave_amplitude.size, int(param.num_timesteps_prop)))
pol_elec = np.zeros((param.microwave_amplitude.size, int(param.num_timesteps_prop)))

# Start timer
start = time.time()

# Calculate system dynamics, looping over microwave amplitude array
for count in range(0, param.microwave_amplitude.size):

    pol_nuc[count], pol_elec[count], pol_i_z_rot, pol_s_z_rot, energies = \
        solid_effect.dynamics(param.microwave_amplitude[count])

    print('{}{:d}{}{:d}{}{:.2f}{}'.format('Finished loop ', (count + 1), ' of ', param.microwave_amplitude.size,
                                          ', total elapsed time ', (time.time() - start), ' s.'))
# End timer
end = time.time()

# Dynamically assign and create output directory
directory = '{}{:.2f}'.format('out/solid_effect/mw_', param.microwave_amplitude[-1]/1E6)
if not os.path.exists(directory):
    os.makedirs(directory)

# Save data and copy parameters file
np.savetxt('{}{}'.format(directory, '/pol_nuc.csv'), pol_nuc, fmt='%.8f', newline='\n')
np.savetxt('{}{}'.format(directory, '/pol_elec.csv'), pol_elec, fmt='%.8f', newline='\n')

# High precision required for sub rotor dynamics (negligible change in file size)
np.savetxt('{}{}'.format(directory, '/pol_i_z_rot.csv'), pol_i_z_rot, fmt='%.12f', newline='\n')
np.savetxt('{}{}'.format(directory, '/pol_s_z_rot.csv'), pol_s_z_rot, fmt='%.12f', newline='\n')
np.savetxt('{}{}'.format(directory, '/energies.csv'), energies, fmt='%.0f', newline='\n')

copyfile("parameters.py", '{}{}'.format(directory, '/parameters.py'))

# Plot data using plotting function
if __name__ == "__main__":
    solid_effect_plotting.plot_all(directory)
    print('End of program.')
    plt.show()
