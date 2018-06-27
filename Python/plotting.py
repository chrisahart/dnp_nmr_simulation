import numpy as np
import matplotlib.pyplot as plt
import parameters as param


def plot_all():

    # Read data from file
    pol_nuc = np.loadtxt('out/pol_nuc.csv')
    pol_elec = np.loadtxt('out/pol_elec.csv')

    # Variables for plotting
    pol_nuc_max = np.zeros(param.microwave_amplitude.size)
    for count in range(0, param.microwave_amplitude.size):
        pol_nuc_max[count] = abs(pol_nuc[count, -1])

    time_array = np.linspace(0, (param.nrot - 1) * (1 / param.freq_rotor), num=param.nrot)
    enhancement_nuc = pol_nuc_max / pol_nuc_max[0]

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
    # fig_pol_nuc_mw.savefig('{}{}{}'.format(filename, '_relax', '.png'), dpi=self.dpi, bbox_inches='tight')

    # Plot enhancement against time
    if pol_nuc.shape[0] >= 20:
        fig_pol_nuc = plt.figure()
        ax_pol_nuc = fig_pol_nuc.add_subplot(111)
        ax_pol_nuc.plot(time_array, abs(pol_nuc[0, :]) / abs(pol_nuc[0, 0]), 'k')
        ax_pol_nuc.plot(time_array, abs(pol_nuc[5, :])/abs(pol_nuc[5, 0]), 'r')
        ax_pol_nuc.plot(time_array, abs(pol_nuc[10, :])/abs(pol_nuc[10, 0]), 'g')
        ax_pol_nuc.plot(time_array, abs(pol_nuc[20, :])/abs(pol_nuc[20, 0]), 'b')
        ax_pol_nuc.set_xlim(time_array[0], time_array[-1])
        ax_pol_nuc.set_ylim(ymin=0)
        ax_pol_nuc.set_xlabel('Time (s)')
        ax_pol_nuc.set_ylabel('Nuclear enhancement')
        fig_pol_nuc.tight_layout()


if __name__ == "__main__":
        plot_all()  # specify directory
        plt.show()
