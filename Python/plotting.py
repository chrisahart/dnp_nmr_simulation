import numpy as np
import matplotlib.pyplot as plt
import parameters as param


def plot_all(directory):

    # Set figure options
    save_dpi = 200

    # Read data from file
    pol_nuc = np.loadtxt('{}{}'.format(directory, '/pol_nuc.csv'))
    pol_elec = np.loadtxt('{}{}'.format(directory, '/pol_elec.csv'))

    # Variables for plotting
    pol_nuc_max = np.zeros(param.microwave_amplitude.size)
    pol_elec_max = np.zeros(param.microwave_amplitude.size)
    for count in range(0, param.microwave_amplitude.size):
        pol_nuc_max[count] = abs(pol_nuc[count, -1])
        pol_elec_max[count] = abs(pol_elec[count, -1])

    time_array = np.linspace(0, (param.nrot - 1) * (1 / param.freq_rotor), num=param.nrot)
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

    # Plot nuclear enhancement against time
    if pol_nuc.shape[0] >= 20:
        fig_pol_nuc = plt.figure()
        ax_pol_nuc = fig_pol_nuc.add_subplot(111)
        ax_pol_nuc.plot(time_array, abs(pol_nuc[0, :]) / abs(pol_nuc[0, 0]), 'k')
        ax_pol_nuc.plot(time_array, abs(pol_nuc[4, :])/abs(pol_nuc[0, 0]), 'r')
        ax_pol_nuc.plot(time_array, abs(pol_nuc[9, :])/abs(pol_nuc[0, 0]), 'g')
        ax_pol_nuc.plot(time_array, abs(pol_nuc[19, :])/abs(pol_nuc[0, 0]), 'b')
        ax_pol_nuc.set_xlim(time_array[0], time_array[-1])
        ax_pol_nuc.set_ylim(ymin=0)
        ax_pol_nuc.set_xlabel('Time (s)')
        ax_pol_nuc.set_ylabel('Nuclear enhancement')
        fig_pol_nuc.tight_layout()
        fig_pol_nuc.savefig('{}{}'.format(directory, '/fig_pol_nuc.png'), dpi=save_dpi, bbox_inches='tight')


if __name__ == "__main__":
        folder = '{}{}'.format('out/solid_effect/microwave_', str(int(param.microwave_amplitude[-1]/1E6)))
        plot_all(folder)
        plt.show()
