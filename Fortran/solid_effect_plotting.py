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
        pol_elec = np.loadtxt('{}{}'.format(directory, '/pol_elec.csv'))

        # Variables for plotting
        pol_nuc_max = np.zeros(param.microwave_amplitude.size)
        pol_elec_max = np.zeros(param.microwave_amplitude.size)
        for count in range(0, param.microwave_amplitude.size):
            pol_nuc_max[count] = abs(pol_nuc[count, -1])
            pol_elec_max[count] = abs(pol_elec[count, -1])

        time_array = np.linspace(0, pol_nuc.shape[1] - 1, num=pol_nuc.shape[1])
        enhancement_nuc = pol_nuc_max / pol_nuc_max[0]
        enhancement_elec = pol_elec_max / pol_elec_max[0]
        enhancement_nuc_time = abs(pol_nuc)

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

        # Plot nuclear enhancement against time
        if param.microwave_amplitude.size >= 20:
            fig_pol_nuc = plt.figure()
            ax_pol_nuc = fig_pol_nuc.add_subplot(111)
            ax_pol_nuc.plot(time_array, enhancement_nuc_time[0, :] / pol_nuc_max[0], 'k',
                            label=param.microwave_amplitude[0]/1E6)
            ax_pol_nuc.plot(time_array, enhancement_nuc_time[4, :] / pol_nuc_max[0], 'r',
                            label=param.microwave_amplitude[4]/1E6)
            ax_pol_nuc.plot(time_array, enhancement_nuc_time[9, :] / pol_nuc_max[0], 'g',
                            label=param.microwave_amplitude[9]/1E6)
            ax_pol_nuc.plot(time_array, enhancement_nuc_time[19, :] / pol_nuc_max[0], 'b',
                            label=param.microwave_amplitude[19]/1E6)
            ax_pol_nuc.legend(loc='upper right', title='Microwave amplitude')
            ax_pol_nuc.set_xlim(time_array[0], time_array[-1])
            ax_pol_nuc.set_ylim(ymin=0)
            ax_pol_nuc.set_xlabel('Time (s)')
            ax_pol_nuc.set_ylabel('Nuclear enhancement')
            fig_pol_nuc.tight_layout()
            fig_pol_nuc.savefig('{}{}'.format(directory, '/nuclear_enhancement_time.png'),
                                dpi=save_dpi, bbox_inches='tight')

    # Plot sub rotor dynamics if data files exist
    elif os.path.exists('{}'.format('out/pol_i_z_rot.out')):

        # Read data from file
        energy_rotor_1 = np.loadtxt("out/energies_1.out") / 1E6
        energy_rotor_2 = np.loadtxt("out/energies_2.out") / 1E6
        energy_rotor_3 = np.loadtxt("out/energies_3.out") / 1E6
        energy_rotor_4 = np.loadtxt("out/energies_4.out") / 1E6
        pol_nuc_rotor = np.loadtxt("out/pol_i_z_rot.out")
        pol_elec_rotor = np.loadtxt("out/pol_s_z_rot.out")
        pol_nuc = np.loadtxt("out/pol_i_z.out")
        pol_elec = np.loadtxt("out/pol_s_z.out")

        # Create plotting arrays
        angle = np.linspace(0, 360, param.time_step_num)
        nuclear_pol = abs(pol_nuc_rotor)/abs(pol_nuc_rotor[0])
        electron_pol = abs(pol_elec_rotor) / abs(pol_elec_rotor[0])
        time_array = np.linspace(0, param.num_timesteps_prop*(1/param.freq_rotor), num=param.num_timesteps_prop)
        enhancement_nuc_time = abs(pol_nuc)
        pol_elec_time = abs(pol_elec)

        # Create 3x3 subplot of polarisation and energy
        label_offset = -0.1
        subfig_fig, (subfig_x1, subfig_x2, subfig_x3) = plt.subplots(3, sharex='all', figsize=(7, 8))

        # Plot electron polarisation
        subfig_x1.plot(angle, electron_pol, 'b')
        subfig_x1.set_ylabel('Elec. pol. / thermal')
        subfig_x1.get_yaxis().set_label_coords(label_offset, 0.5)

        # Plot energy
        subfig_x2.plot(angle, energy_rotor_1, 'k')
        subfig_x2.plot(angle, energy_rotor_2, 'r')
        subfig_x2.plot(angle, energy_rotor_3, 'b')
        subfig_x2.plot(angle, energy_rotor_4, 'g')
        subfig_x2.set_ylabel('Energy / MHz')
        subfig_x2.get_yaxis().set_label_coords(label_offset, 0.5)

        # Plot nuclear polarisation
        subfig_x3.plot(angle, nuclear_pol, 'r')
        subfig_x3.set_ylabel('Nuc. pol. / thermal')
        subfig_x3.get_yaxis().set_label_coords(label_offset, 0.5)

        # Adjust subplot x shared label
        subfig_x3.set_xlim(angle[0], angle[-1])
        subfig_x3.set_xlabel('Rotor position (degrees)')
        subfig_fig.tight_layout()
        subfig_fig.subplots_adjust(hspace=0)

        subfig_fig.savefig('{}{}'.format(directory, '/rotor_dynamics.png'), dpi=save_dpi, bbox_inches='tight')

        # Plot nuclear enhancement against time
        fig_pol_nuc = plt.figure()
        ax_pol_nuc = fig_pol_nuc.add_subplot(111)
        ax_pol_nuc.plot(time_array, enhancement_nuc_time / enhancement_nuc_time[0], 'k')
        ax_pol_nuc.set_xlim(time_array[0], time_array[-1])
        ax_pol_nuc.set_ylim(ymin=0)
        ax_pol_nuc.set_xlabel('Time (s)')
        ax_pol_nuc.set_ylabel('Nuclear enhancement')
        fig_pol_nuc.tight_layout()
        fig_pol_nuc.savefig('{}{}'.format(directory, '/nuclear_enhancement.png'), dpi=save_dpi, bbox_inches='tight')

        # Plot electron polarisation against time
        fig_pol_elec = plt.figure()
        ax_pol_elec = fig_pol_elec.add_subplot(111)
        ax_pol_elec.plot(time_array, pol_elec_time / pol_elec_time[0], 'k')
        ax_pol_elec.set_xlim(time_array[0], time_array[-1])
        ax_pol_elec.set_ylim(ymin=0)
        ax_pol_elec.set_xlabel('Time (s)')
        ax_pol_elec.set_ylabel('Electron polarisation')
        fig_pol_elec.tight_layout()
        fig_pol_elec.savefig('{}{}'.format(directory, '/electron_polarisation.png'), dpi=save_dpi, bbox_inches='tight')

    print('Finished plotting.')


if __name__ == "__main__":
        folder = 'out'
        plot_all(folder)
        plt.show()