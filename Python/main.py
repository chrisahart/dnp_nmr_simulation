import numpy as np
import matplotlib.pyplot as plt
import dynamics as dynamics
import time
import parameters as param

# TODO Finish adding dynamics code into functions, add output to file and plotting code

""" Call dynamics function as required, saving output to working directory.
    Powder averaging, and effect of microwave amplitude can be simulated.
"""

# Input parameters
microwave_amplitude = np.arange(1, 20, 5) * 1E6
nrot = int(np.round(40 * param.freq_rotor))  # TODO Figure out what this does

# Pre-allocate arrays
pol_nuc = np.zeros((microwave_amplitude.size, int(nrot)))
pol_nuc_max = np.zeros(microwave_amplitude.size)
pol_elec = np.zeros((microwave_amplitude.size, int(nrot)))
pol_elec_max = np.zeros(microwave_amplitude.size)

start = time.time()

for count in range(0, microwave_amplitude.size):

    pol_nuc[count], pol_elec[count] = dynamics.dynamics(microwave_amplitude[count])
    pol_nuc_max[count] = abs(pol_nuc[count, -1])
    pol_elec_max[count] = abs(pol_elec[count, -1])

    print('Finished loop', (count + 1), 'of', microwave_amplitude.size, '.')  # TODO improve
    print('Elapsed time:', np.round(time.time() - start), ' s. \n')

end = time.time()

# Variables for plotting
time_array = np.linspace(0, (nrot-1)*(1/param.freq_rotor), num=nrot)
enhancement_nuc = pol_nuc_max / pol_nuc_max[0]
# TODO create loop to calculate enhancement nuc for plotting against time

# Plot enhancement against microwave amplitude TODO make prettier
fig_pol_nuc_mw = plt.figure()
ax_pol_nuc_mw = fig_pol_nuc_mw.add_subplot(111)
ax_pol_nuc_mw.plot(microwave_amplitude/1E6, enhancement_nuc,
                   'kx', microwave_amplitude/1E6, enhancement_nuc, 'k')
ax_pol_nuc_mw.set_xlim(microwave_amplitude[0]/1E6, microwave_amplitude[-1]/1E6)
ax_pol_nuc_mw.set_ylim(left=0)
ax_pol_nuc_mw.set_xlabel('Microwave amplitude / MHz')
ax_pol_nuc_mw.set_ylabel('Nuclear enhancement')
fig_pol_nuc_mw.tight_layout()

# Plot enhancement against time
fig_pol_nuc = plt.figure()
ax_pol_nuc = fig_pol_nuc.add_subplot(111)
ax_pol_nuc.plot(time_array, abs(pol_nuc[0, :])/abs(pol_nuc[0, 0]), 'k')
ax_pol_nuc.plot(time_array, abs(pol_nuc[5, :])/abs(pol_nuc[5, 0]), 'r')
ax_pol_nuc.plot(time_array, abs(pol_nuc[10, :])/abs(pol_nuc[10, 0]), 'g')
ax_pol_nuc.plot(time_array, abs(pol_nuc[20, :])/abs(pol_nuc[20, 0]), 'b')
ax_pol_nuc.set_xlim(time_array[0], time_array[-1])
ax_pol_nuc.set_ylim(left=0)
ax_pol_nuc.set_xlabel('Time (s)')
ax_pol_nuc.set_ylabel('Nuclear enhancement')
fig_pol_nuc.tight_layout()

plt.show()
