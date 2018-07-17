import numpy as np
import matplotlib.pyplot as plt
import parameters as param
import os


def plot_all(directory):

    # Set figure options
    save_dpi = 200

    # Plot multi rotor dynamics if array size greater than one
    if param.microwave_amplitude.size > 1:

        # Read data from file
        pol_nuc = np.loadtxt('{}{}'.format(directory, '/pol_nuc.csv'))
        pol_elec1 = np.loadtxt('{}{}'.format(directory, '/pol_elec1.csv'))
        pol_elec2 = np.loadtxt('{}{}'.format(directory, '/pol_elec2.csv'))

        # Variables for plotting
        pol_nuc_max = np.zeros(param.microwave_amplitude.size)
        pol_elec_max = np.zeros(param.microwave_amplitude.size)
        for count in range(0, param.microwave_amplitude.size):
            pol_nuc_max[count] = abs(pol_nuc[count, -1])
            pol_elec_max[count] = abs(pol_elec1[count, -1]) + abs(pol_elec2[count, -1])

        enhancement_nuc = pol_nuc_max / pol_nuc_max[0]
        enhancement_elec = pol_elec_max / pol_elec_max[0]

        # Plot nuclear enhancement against microwave amplitude
        fig_pol_nuc_mw = plt.figure()
        ax_pol_nuc_mw = fig_pol_nuc_mw.add_subplot(111)
        ax_pol_nuc_mw.plot(param.microwave_amplitude / 1E6, enhancement_nuc,
                           'kx', param.microwave_amplitude / 1E6, enhancement_nuc, 'k')
        ax_pol_nuc_mw.set_xlim(0, param.microwave_amplitude[-1] / 1E6)
        ax_pol_nuc_mw.set_ylim(ymin=0)
        ax_pol_nuc_mw.set_xlabel('Microwave amplitude / MHz')
        ax_pol_nuc_mw.set_ylabel('Nuclear enhancement')
        fig_pol_nuc_mw.tight_layout()
        fig_pol_nuc_mw.savefig('{}{}'.format(directory, '/fig_pol_nuc_mw.png'), dpi=save_dpi, bbox_inches='tight')

        # Plot electron enhancement against microwave amplitude
        fig_pol_elec_mw = plt.figure()
        ax_pol_elec_mw = fig_pol_elec_mw.add_subplot(111)
        ax_pol_elec_mw.plot(param.microwave_amplitude / 1E6, enhancement_elec,
                            'kx', param.microwave_amplitude / 1E6, enhancement_elec, 'k')
        ax_pol_elec_mw.set_xlim(0, param.microwave_amplitude[-1] / 1E6)
        ax_pol_elec_mw.set_ylim(ymin=0)
        ax_pol_elec_mw.set_xlabel('Microwave amplitude / MHz')
        ax_pol_elec_mw.set_ylabel('Electron enhancement')
        fig_pol_elec_mw.tight_layout()
        fig_pol_elec_mw.savefig('{}{}'.format(directory, '/fig_pol_elec_mw.png'), dpi=save_dpi, bbox_inches='tight')

    # Plot sub rotor dynamics if data files exist
    elif os.path.exists('{}{}'.format(directory, '/pol_i_z_rot.csv')):

        # Read data from file
        pol_nuc_rotor = np.loadtxt('{}{}'.format(directory, '/pol_i_z_rot.csv'))
        pol_elec_rotor1 = np.loadtxt('{}{}'.format(directory, '/pol_s1_z_rot.csv'))
        pol_elec_rotor2 = np.loadtxt('{}{}'.format(directory, '/pol_s2_z_rot.csv'))
        energy_rotor = np.loadtxt('{}{}'.format(directory, '/energies.csv'))
        pol_nuc = np.loadtxt('{}{}'.format(directory, '/pol_nuc.csv'))
        pol_elec1 = np.loadtxt('{}{}'.format(directory, '/pol_elec1.csv'))
        pol_elec2 = np.loadtxt('{}{}'.format(directory, '/pol_elec1.csv'))

        # Create plotting arrays
        angle = np.linspace(0, 360, param.time_step_num)
        nuclear_pol = abs(pol_nuc_rotor)/abs(pol_nuc_rotor[0])
        electron_pol1 = abs(pol_elec_rotor1) / abs(pol_elec_rotor1[0])
        electron_pol2 = abs(pol_elec_rotor2) / abs(pol_elec_rotor2[0])
        energy_rotor = energy_rotor / 1E6
        time_array = np.linspace(0, param.num_periods, num=param.num_timesteps_prop)
        enhancement_nuc_time = abs(pol_nuc)
        pol_elec1_time = abs(pol_elec1)
        pol_elec2_time = abs(pol_elec2)
        pol_elec_time = pol_elec1_time + pol_elec2_time

        # Create 3x3 subplot of polarisation and energy
        label_offset = -0.1
        subfig_fig, (subfig_x1, subfig_x2, subfig_x3, subfig_x4) = plt.subplots(4, sharex='all', figsize=(7, 9))

        # Plot electron polarisation
        subfig_x1.plot(angle, electron_pol1, 'b')
        subfig_x1.set_ylabel('Elec. 1 pol. / thermal')
        subfig_x1.get_yaxis().set_label_coords(label_offset, 0.5)

        # Plot electron polarisation
        subfig_x2.plot(angle, electron_pol2, 'b')
        subfig_x2.set_ylabel('Elec. 2 pol. / thermal')
        subfig_x2.get_yaxis().set_label_coords(label_offset, 0.5)

        # Plot energy
        subfig_x3.plot(angle, energy_rotor[:, 0], 'k')
        subfig_x3.plot(angle, energy_rotor[:, 1], 'r')
        subfig_x3.plot(angle, energy_rotor[:, 2], 'b')
        subfig_x3.plot(angle, energy_rotor[:, 3], 'g')
        subfig_x3.plot(angle, energy_rotor[:, 4], 'orange')
        subfig_x3.plot(angle, energy_rotor[:, 5], 'm')
        subfig_x3.plot(angle, energy_rotor[:, 6], 'brown')
        subfig_x3.plot(angle, energy_rotor[:, 7], 'grey')
        subfig_x3.set_ylabel('Energy / MHz')
        subfig_x3.get_yaxis().set_label_coords(label_offset, 0.5)

        # Plot nuclear polarisation
        subfig_x4.plot(angle, nuclear_pol, 'r')
        subfig_x4.set_ylabel('Nuc. pol. / thermal')
        subfig_x4.get_yaxis().set_label_coords(label_offset, 0.5)

        # Adjust subplot x shared label
        subfig_x4.set_xlim(angle[0], angle[-1])
        subfig_x4.set_xlabel('Rotor position (degrees)')
        subfig_fig.tight_layout()
        subfig_fig.subplots_adjust(hspace=0)

        subfig_fig.savefig('{}{}'.format(directory, '/rotor_dynamics.png'), dpi=save_dpi, bbox_inches='tight')

        # Create 2x2 subplot of multi-rotor polarisation
        subfig2_fig, (subfig2_x1, subfig2_x2) = plt.subplots(2, sharex='all', figsize=(6, 8))

        # Plot electron polarisation
        subfig2_x1.plot(time_array, pol_elec_time / pol_elec_time[0], 'b')
        subfig2_x1.set_ylim(0, 1)
        subfig2_x1.set_ylabel('Elec. pol. / thermal')

        # Plot nuclear polarisation
        subfig2_x2.plot(time_array, enhancement_nuc_time / enhancement_nuc_time[0], 'r')
        subfig2_x2.set_ylim(ymin=0)
        subfig2_x2.set_ylabel('Nuc. pol. / thermal')

        # Adjust subplot x shared label
        subfig2_x2.set_xlim(time_array[0] - 0.1, time_array[-1])
        subfig2_x2.set_xlabel('Time (s)')
        subfig2_fig.tight_layout()
        subfig2_fig.subplots_adjust(hspace=0)

        subfig2_fig.savefig('{}{}'.format(directory, '/multi_rotor.png'), dpi=save_dpi, bbox_inches='tight')

    print('Finished plotting.')


if __name__ == "__main__":
        folder = '{}{:.2f}'.format('out/cross_effect/mw_', param.microwave_amplitude[-1] / 1E6)
        plot_all(folder)
        plt.show()
