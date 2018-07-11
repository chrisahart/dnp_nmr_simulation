program solid_effect_main

    use omp_lib
    use f2py_dynamics
    use f2py_functions
    implicit none

    real(kind = 8), parameter :: PI = 4.D0 * DATAN(1.D0)
    real(kind = 8), parameter :: deg_to_rad = PI / 180.D0

    ! todo make these allocatable based on time_num
    real(kind = 8), dimension(120000)  :: pol_i_z, pol_s_z
    real(kind = 8), dimension(10000) :: pol_i_z_rot, pol_s_z_rot
    real(kind = 8), dimension(10000, 4) :: energies

    integer :: time_num, time_num_prop, num_periods
    real(kind = 8), dimension(3) :: gtensor, hyperfine_angles, orientation_se
    real(kind = 8) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
    real(kind = 8) :: time_step, freq_rotor, b_field, microwave_amplitude, temperature
    real(kind = 8) :: t1_nuc, t1_elec, t2_nuc, t2_elec
    real(kind = 8) :: wtime

    ! Literature parameters todo match python parameters file exactly
    b_field = 9.4D0                                                       ! Magnetic field
    temperature = 100D0                                                   ! Temperature
    nuclear_frequency = -400.9D6                                           ! Nuclear frequency 1
    electron_frequency = 28.025D9 * b_field                             ! Electron frequency
    microwave_frequency = 264D9                                         ! Microwave frequency
    freq_rotor = 3D3                                                    ! Rotor frequency
    orientation_se = (/253.6D0, 105.1D0, 123.8D0   /) * deg_to_rad !G anisotropy angles for electron 1 (SE)
    gtensor = (/(2.00614D0 / 2.D0), (2.00194D0 / 2.D0), (2.00988D0 / 2.D0) /) ! G tensor principal values
    hyperfine_coupling = 3D6                                            ! Hyperfine coupling amplitude
    hyperfine_angles = (/0.D0, 0.D0, 0.D0 /)                         ! Hyperfine angles for e1-n
    t1_elec = 0.3D-3                                                    ! Electron spin-lattice relaxation T1 (s)
    t1_nuc = 10.D0                                                         ! Nuclear spin-lattice relaxation T1 (s)
    t2_elec = 1D-6                                                      ! T2 electron
    t2_nuc = 1D-3                                                       ! T2 nucleus
    microwave_amplitude = 0.85D0 * 1D6                        ! Microwave field amplitude np.arange(1, 10, 5)

    ! System variables
    time_num = 1D4                                            ! Number of timesteps within rotor period
    time_step = (1.D0 / freq_rotor) / time_num                        ! Value of timestep within rotor period
    num_periods = 40                                                    ! Number of rotor periods to propagate system
    time_num_prop = num_periods * int(freq_rotor)      ! Number of timesteps to propagate system

    wtime = omp_get_wtime()

    call main(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, &
            orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, microwave_amplitude, &
            t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, time_num_prop, &
            pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies)

    wtime = omp_get_wtime () - wtime
    write(6, *) sngl(wtime)

    ! Save output, todo replace with FLIBS to save to csv, and dynamic save location
    open(1, file = 'out/pol_i_z.out', status = 'replace')
    open(2, file = 'out/pol_s_z.out', status = 'replace')
    open(3, file = 'out/pol_i_z_rot.out', status = 'replace')
    open(4, file = 'out/pol_s_z_rot.out', status = 'replace')
    open(5, file = 'out/energies_1.out', status = 'replace')
    open(7, file = 'out/energies_2.out', status = 'replace')
    open(8, file = 'out/energies_3.out', status = 'replace')
    open(9, file = 'out/energies_4.out', status = 'replace')

    write(1, *) pol_i_z
    write(2, *) pol_s_z
    write(3, *) pol_i_z_rot
    write(4, *) pol_s_z_rot
    write(5, *) energies(:, 1)
    write(7, *) energies(:, 2)
    write(8, *) energies(:, 3)
    write(9, *) energies(:, 4)

    close(1, status = 'keep')
    close(2, status = 'keep')
    close(3, status = 'keep')
    close(4, status = 'keep')
    close(5, status = 'keep')
    close(7, status = 'keep')
    close(8, status = 'keep')
    close(9, status = 'keep')

    write(6,*) 'Finished saving data, end of program.'

end program solid_effect_main