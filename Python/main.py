import numpy as np
import matplotlib.pyplot as plt
import dynamics as dynamics
import time

# TODO Finish adding dynamics code into functions, add output to file and plotting code

""" Powder averaging, call dynamics function as required, saving output to working directory
"""

start = time.time()

energies, iznew, sznew = dynamics.dynamics()

end = time.time()
print('Elapsed time', end - start)

fig_energy = plt.figure()
ax_energy = fig_energy.add_subplot(111)
ax_energy.plot(energies[:, 0])
ax_energy.plot(energies[:, 1])
ax_energy.plot(energies[:, 2])
ax_energy.plot(energies[:, 3])

fig_pol_nuc = plt.figure()
ax_pol_nuc = fig_pol_nuc.add_subplot(111)
ax_pol_nuc.plot(abs(iznew))

fig_pol_elec = plt.figure()
ax_pol_elec = fig_pol_elec.add_subplot(111)
ax_pol_elec.plot(abs(sznew))
ax_pol_elec.set_xlim(0, 20)

plt.show()
