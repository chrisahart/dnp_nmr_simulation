program solid_effect_main

    use omp_lib
    use f2py_dynamics
    use f2py_functions
    implicit none

    integer, intent(in) :: time_num, time_num_prop
    real(kind = 8), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_se
    real(kind = 8), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
    real(kind = 8), intent(in) :: time_step, freq_rotor
    real(kind = 8), intent(in) :: microwave_amplitude
    real(kind = 8), intent(in) :: t1_nuc, t1_elec, t2_nuc, t2_elec, temperature

    real(kind = 8), dimension(time_num_prop), intent(out) :: pol_i_z, pol_s_z
    real(kind = 8), dimension(time_num), intent(out) :: pol_i_z_rot, pol_s_z_rot
    real(kind = 8), dimension(time_num, 4), intent(out) :: energies

    double precision wtime

    wtime = omp_get_wtime()

    call main(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, &
            orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, microwave_amplitude, &
            t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, time_num_prop, &
            pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies)

    wtime = omp_get_wtime () - wtime
    write(6, *) sngl(wtime)

    ! Save output

end program solid_effect_main