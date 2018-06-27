import numpy as np
import matplotlib.pyplot as plt
import dynamics as dynamics
import time
import parameters as param

# TODO Finish adding dynamics code into functions, add output to file and plotting code

""" Call dynamics function as required, saving output to working directory.
    Powder averaging, and effect of microwave amplitude can be simulated.
"""

microwave_amplitude = np.arange(1, 200, 5) * 1E6

nrot = int(np.round(40 * param.freq_rotor))
iznew = np.zeros((microwave_amplitude.size, int(nrot)))  # TODO generalize
iznew_max = np.zeros(microwave_amplitude.size)  # TODO generalize
sznew = np.zeros((microwave_amplitude.size, int(nrot)))  # TODO generalize
sznew_max = np.zeros(microwave_amplitude.size)  # TODO generalize
energies = np.zeros((int(param.time_step_num), 4), dtype=np.complex)  # TODO generalize

start = time.time()

for count in range(0, microwave_amplitude.size):

    iznew[count], sznew[count] = dynamics.dynamics(microwave_amplitude[count])
    iznew_max[count] = abs(iznew[count, -1])
    sznew_max[count] = abs(sznew[count, -1])
    print('Finished loop', (count + 1), 'of', microwave_amplitude.size)

end = time.time()
print('Elapsed time', end - start)

time = np.linspace(0, nrot-1, num=nrot)

# Plot enhancement against microwave amplitude
fig_pol_nuc_mw = plt.figure()
ax_pol_nuc_mw = fig_pol_nuc_mw.add_subplot(111)
ax_pol_nuc_mw.plot(microwave_amplitude, iznew_max, 'kx', microwave_amplitude, iznew_max, 'k')
ax_pol_nuc_mw.set_xlabel('microwave_amplitude?')
ax_pol_nuc_mw.set_ylabel('enhancement')
fig_pol_nuc_mw.tight_layout()

# Plot enhancement against time
fig_pol_nuc = plt.figure()
ax_pol_nuc = fig_pol_nuc.add_subplot(111)
ax_pol_nuc.plot(time, abs(iznew[0, :]), 'r')
ax_pol_nuc.plot(time, abs(iznew[1, :]), 'g')
ax_pol_nuc.plot(time, abs(iznew[2, :]), 'b')
ax_pol_nuc.plot(time, abs(iznew[3, :]), 'y')
ax_pol_nuc.set_xlabel('time?')
ax_pol_nuc.set_ylabel('enhancement')
fig_pol_nuc.tight_layout()

# Plot energy against time
# fig_energy = plt.figure()
# ax_energy = fig_energy.add_subplot(111)
# ax_energy.plot(energies[:, 0])
# ax_energy.plot(energies[:, 1])
# ax_energy.plot(energies[:, 2])
# ax_energy.plot(energies[:, 3])

plt.show()
