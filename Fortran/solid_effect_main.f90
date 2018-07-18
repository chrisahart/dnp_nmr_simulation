program solid_effect_main

    ! Call dynamics module to calculate dynamics for solid effect MAS DNP NMR.

    use omp_lib
    use solid_effect_dynamics
    use iso_fortran_env
    implicit none

    integer, parameter :: wp = real64
    integer, parameter :: r15 = selected_real_kind(15)
    real(kind = 8), parameter :: PI = 4.D0 * DATAN(1.D0)
    real(kind = 8), parameter :: rad = PI / 180.D0
    integer, parameter :: sizeH = 2 ** (2)

    real(kind = 8), allocatable :: energies(:, :), pol_i_z(:), pol_s_z(:), pol_i_z_rot(:), pol_s_z_rot(:)
    integer :: time_num, time_num_prop, num_periods
    real(kind = 8), dimension(3) :: gtensor, hyperfine_angles, orientation_se, orientation_ce_1, orientation_ce_2
    real(kind = 8) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
    real(kind = 8) :: time_step, freq_rotor, b_field, microwave_amplitude, temperature
    real(kind = 8) :: t1_nuc, t1_elec, t2_nuc, t2_elec
    real(kind = 8) :: wtime

    ! Literature parameters
    b_field = 9.4D0                                                           ! Magnetic field
    temperature = 100D0                                                       ! Temperature
    nuclear_frequency = -400.9D6                                              ! Nuclear frequency 1
    electron_frequency = 28.025D9 * b_field                                   ! Electron frequency
    microwave_frequency = 264D9                                               ! Microwave frequency
    freq_rotor = 3D3                                                          ! Rotor frequency
    orientation_se = [253.6D0, 105.1D0, 123.8D0] * rad                        ! G anisotropy angles for electron 1 (SE)
    orientation_ce_1 = [253.6D0, 105.1D0, 123.8D0] * rad                      ! G anisotropy angles for electron 1 (CE)
    orientation_ce_2 = orientation_ce_1 + [102.D0, 104.D0, 124.D0] * rad      ! G anisotropy angles for electron 2 (CE)
    gtensor = [(2.00614D0 / 2.D0), (2.00194D0 / 2.D0), (2.00988D0 / 2.D0)]    ! G tensor principal values
    hyperfine_coupling = 3D6                                                  ! Hyperfine coupling amplitude
    hyperfine_angles = [0.D0, 0.D0, 0.D0]                                     ! Hyperfine angles for e1-n
    t1_elec = 0.3D-3                                                          ! Electron spin-lattice relaxation T1 (s)
    t1_nuc = 10.D0                                                            ! Nuclear spin-lattice relaxation T1 (s)
    t2_elec = 1D-6                                                            ! T2 electron
    t2_nuc = 1D-3                                                             ! T2 nucleus
    microwave_amplitude = 0.85D0 * 1D6                                        ! Microwave field amplitude

    ! System variables
    time_num = 1D4                                                            ! Number of timesteps within rotor period
    time_step = (1.D0 / freq_rotor) / time_num                                ! Value of timestep within rotor period
    num_periods = 40                                                          ! Number of rotor periods
    time_num_prop = num_periods * int(freq_rotor)                             ! Number of timesteps to propagate system

    ! Allocate output variables based on input parameters
    allocate (energies(time_num, sizeH), pol_i_z(time_num_prop), pol_s_z(time_num_prop))
    allocate (pol_i_z_rot(time_num), pol_s_z_rot(time_num))

    ! Manually set number of OMP threads
    call omp_set_num_threads(8)

    ! Start timer (using OpenMP to work across multiple cores)
    wtime = omp_get_wtime()

    ! Call main() to calculate SE dynamics
    call main(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, &
            orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, microwave_amplitude, &
            t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, time_num_prop, &
            pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies)

    wtime = omp_get_wtime () - wtime
    write(6, *) 'Total elapsed time:', sngl(wtime)

    ! Open output files
    open(1, file = 'out/pol_i_z.out', status = 'replace')
    open(2, file = 'out/pol_s_z.out', status = 'replace')
    open(3, file = 'out/pol_i_z_rot.out', status = 'replace')
    open(4, file = 'out/pol_s_z_rot.out', status = 'replace')
    open(5, file = 'out/energies_1.out', status = 'replace')
    open(7, file = 'out/energies_2.out', status = 'replace')
    open(8, file = 'out/energies_3.out', status = 'replace')
    open(9, file = 'out/energies_4.out', status = 'replace')

    ! Write to output files
    write(1, *) pol_i_z
    write(2, *) pol_s_z
    write(3, *) pol_i_z_rot
    write(4, *) pol_s_z_rot
    write(5, *) energies(:, 1)
    write(7, *) energies(:, 2)
    write(8, *) energies(:, 3)
    write(9, *) energies(:, 4)

    ! Close output files
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