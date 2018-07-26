program solid_effect_main

    ! Call dynamics module to calculate dynamics for solid effect MAS DNP NMR.

    use omp_lib
    use solid_effect_dynamics
    use iso_fortran_env
    implicit none

    integer, parameter :: wp = selected_real_kind(15, 307)
    real(wp), parameter :: pi = 4._wp * atan(1._wp)
    real(wp), parameter :: rad = pi / 180._wp
    integer, parameter :: sizeH = 2 ** (2)

    real(wp), allocatable :: energies(:, :), pol_i_z(:), pol_s_z(:), pol_i_z_rot(:), pol_s_z_rot(:), pol_iz_final(:)
    real(wp), dimension(3) :: gtensor, hyperfine_angles, orientation_se, orientation_ce_1, orientation_ce_2
    real(wp) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
    real(wp) :: time_step, freq_rotor, b_field, temperature
    real(wp) :: t1_nuc, t1_elec, t2_nuc, t2_elec, wtime
    integer :: time_num, time_num_prop, num_periods, count

    ! Literature parameters
    real(wp) :: microwave_amplitude(40) = [(count, count = 1, 200, 5)] * 1E6          ! Microwave amplitude array
    !real(wp) :: microwave_amplitude(1) = [0.85E6_wp]                             ! Microwave amplitude
    b_field = 9.4_wp                                                           ! Magnetic field
    temperature = 100_wp                                                       ! Temperature
    nuclear_frequency = -400.9E6_wp                                               ! Nuclear frequency 1
    electron_frequency = 28.025E9_wp * b_field                                    ! Electron frequency
    microwave_frequency = 264E9_wp !265.2E9_wp !  ! !  !                                        ! Microwave frequency
    freq_rotor = 3E3_wp                                                           ! Rotor frequency
    orientation_se = [253.6_wp, 105.1_wp, 123.8_wp] * rad                      ! G anisotropy angles for electron 1 (SE)
    orientation_ce_1 = [253.6_wp, 105.1_wp, 123.8_wp] * rad                    ! G anisotropy angles for electron 1 (CE)
    orientation_ce_2 = orientation_ce_1 + [102._wp, 104._wp, 124._wp] * rad    ! G anisotropy angles for electron 2 (CE)
    gtensor = [(2.00614_wp / 2._wp), (2.00194_wp / 2._wp), (2.00988_wp / 2._wp)]        ! G tensor principal values
    hyperfine_coupling = 3E6_wp                                                   ! Hyperfine coupling amplitude
    hyperfine_angles = [0._wp, 0._wp, 0._wp] ![253.6, 105.1, 123.8] * rad !Hyperfine angles for e1-n
    t1_elec = 0.3E-3_wp                                                        ! Electron spin-lattice relaxation T1 (s)
    t1_nuc = 10._wp                                                            ! Nuclear spin-lattice relaxation T1 (s)
    t2_elec = 1E-6_wp                                                             ! T2 electron
    t2_nuc = 1E-3_wp                                                              ! T2 nucleus

    ! System variables
    time_num = 1E4                                                             ! Number of timesteps within rotor period
    time_step = (1._wp / freq_rotor) / time_num                                ! Value of timestep within rotor period
    num_periods = 40                                                           ! Number of rotor periods
    time_num_prop = num_periods * int(freq_rotor)                              ! Number of timesteps to propagate system

    ! Allocate output variables based on input parameters
    allocate (energies(time_num, sizeH), pol_i_z(time_num_prop), pol_s_z(time_num_prop))
    allocate (pol_i_z_rot(time_num), pol_s_z_rot(time_num), pol_iz_final(size(microwave_amplitude)))

    ! Manually set number of OMP threads
    call omp_set_num_threads(8)

    ! Start timer (using OpenMP to work across multiple cores)
    wtime = omp_get_wtime()

    !!$omp parallel do default(private) &
    !!$omp& shared(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, orientation_se) &
    !!$omp& shared(electron_frequency, microwave_frequency, nuclear_frequency, microwave_amplitude, t1_nuc, t1_elec) &
    !!$omp& shared(t2_nuc, t2_elec, temperature, time_num_prop, pol_iz_final)
    do count = 1, size(microwave_amplitude)

        ! Call main() to calculate SE dynamics
        call main(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, &
                orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, &
                 microwave_amplitude(count), t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, time_num_prop, &
                pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies)

        pol_iz_final(count) = abs(pol_i_z(time_num_prop))
        write(6, *) 'Finished loop', count, 'of',  size(microwave_amplitude), '.'

    end do
    !!$omp end parallel do

    wtime = omp_get_wtime () - wtime
    write(6, *) 'Total elapsed time:', sngl(wtime)

    ! Open output files
    open(1, file = 'out/pol_i_z.out', status = 'replace')
    open(2, file = 'out/pol_s_z.out', status = 'replace')
    open(3, file = 'out/pol_i_z_rot.out', status = 'replace')
    open(4, file = 'out/pol_s_z_rot.out', status = 'replace')
    open(5, file = 'out/microwave_amplitude.out', status = 'replace')
    open(7, file = 'out/num_periods.out', status = 'replace')
    open(8, file = 'out/pol_iz_final.out', status = 'replace')
    open(9, file = 'out/energies_1.out', status = 'replace')
    open(10, file = 'out/energies_2.out', status = 'replace')
    open(11, file = 'out/energies_3.out', status = 'replace')
    open(12, file = 'out/energies_4.out', status = 'replace')

    ! Write to output files
    write(1, *) pol_i_z
    write(2, *) pol_s_z
    write(3, *) pol_i_z_rot
    write(4, *) pol_s_z_rot
    write(5, *) microwave_amplitude
    write(7, *) num_periods
    write(8, *) pol_iz_final
    write(9, *) energies(:, 1)
    write(10, *) energies(:, 2)
    write(11, *) energies(:, 3)
    write(12, *) energies(:, 4)

    ! Close output files
    close(1, status = 'keep')
    close(2, status = 'keep')
    close(3, status = 'keep')
    close(4, status = 'keep')
    close(5, status = 'keep')
    close(7, status = 'keep')
    close(8, status = 'keep')
    close(9, status = 'keep')
    close(10, status = 'keep')
    close(11, status = 'keep')
    close(12, status = 'keep')

    write(6,*) 'Finished saving data, end of program.'

end program solid_effect_main