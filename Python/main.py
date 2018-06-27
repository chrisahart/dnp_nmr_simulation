import numpy as np
import matplotlib.pyplot as plt
import dynamics as dynamics
import time
import parameters as param
import plotting
from shutil import copyfile
import os

""" Call dynamics function as required, saving output to disk.
    Powder averaging, and effect of microwave amplitude can be simulated.
"""

# Pre-allocate arrays
pol_nuc = np.zeros((param.microwave_amplitude.size, int(param.nrot)))
pol_elec = np.zeros((param.microwave_amplitude.size, int(param.nrot)))

# Start timer
start = time.time()

# Calculate system dynamics, looping over microwave amplitude array
for count in range(0, param.microwave_amplitude.size):

    pol_nuc[count], pol_elec[count] = dynamics.dynamics(param.microwave_amplitude[count])

    print('{}{:d}{}{:d}{}{:.0f}{}'.format('Finished loop ', (count + 1), ' of ', param.microwave_amplitude.size,
                                          ', total elapsed time ', np.round(time.time() - start), ' s.'))
# End timer
end = time.time()

# Dynamically assign and create output directory
directory = '{}{}'.format('out/solid_effect/microwave_', str(int(param.microwave_amplitude[-1]/1E6)))
if not os.path.exists(directory):
    os.makedirs(directory)

# Save data and copy parameters file
np.savetxt('{}{}'.format(directory, '/pol_nuc.csv'), pol_nuc, fmt='%.8f', newline='\n')
np.savetxt('{}{}'.format(directory, '/pol_elec.csv'), pol_elec, fmt='%.8f', newline='\n')
copyfile("parameters.py", '{}{}'.format(directory, '/parameters.py'))

# Plot data using plotting function
if __name__ == "__main__":
    plotting.plot_all(directory)
    plt.show()
