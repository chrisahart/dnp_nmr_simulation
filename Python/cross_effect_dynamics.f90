module cross_effect_dynamics

    ! This module calculates dynamics for solid effect MAS DNP NMR.

contains

    subroutine main(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, orientation_ce_1, &
                orientation_ce_2, electron_frequency, microwave_frequency, nuclear_frequency, microwave_amplitude, &
                t1_nuc,  t1_elec, t2_nuc, t2_elec, temperature, time_num_prop, &
                pol_i_z, pol_s1_z, pol_s2_z, pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot, energies)

        ! Call required functions to calculate dynamics for solid effect MAS DNP NMR

        use omp_lib
        use functions
        use iso_fortran_env
        implicit none

        integer, parameter :: sizeH = 2 ** (3), sizeL = 4 ** (3) ! Change (value) to change number of spins
        integer, parameter :: wp = selected_real_kind(15, 307)  ! Standard double precision

        integer, intent(in):: time_num, time_num_prop
        real(wp), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_ce_1, orientation_ce_2
        real(wp), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(wp), intent(in) :: time_step, freq_rotor, microwave_amplitude
        real(wp), intent(in) :: t1_nuc, t1_elec, t2_nuc, t2_elec, temperature

        real(wp), dimension(time_num_prop), intent(out) :: pol_i_z, pol_s1_z, pol_s2_z
        real(wp), dimension(time_num), intent(out) :: pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot
        real(wp), dimension(time_num, sizeH), intent(out) :: energies

        complex(wp), dimension(time_num, sizeL, sizeL)  :: propagator
        real(wp), dimension(sizeH, sizeH) :: density_mat
        real(wp) :: eig_vector(time_num, sizeH, sizeH), eig_vector_inv(time_num, sizeH, sizeH)
        real(wp) :: eigvals(time_num, sizeH), eigvectors_temp(time_num, sizeH, sizeH)
        real(wp), dimension(time_num, sizeH, sizeH) :: hamiltonian
        integer :: count1, count2
        integer(wp) :: indices(sizeH)

        !real(wp) :: wtime

        ! Construct intrinsic Hilbert space Hamiltonian
        call calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, temperature, hyperfine_coupling, &
                hyperfine_angles, orientation_ce_1, orientation_ce_2, electron_frequency, microwave_frequency, &
                nuclear_frequency, sizeH, sizeL, hamiltonian, density_mat)

        ! Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
        call eig_real(hamiltonian, eigvals, eig_vector)

        ! Sort eigenvalues and eigenvectors at zero time (only adds around 1ms to total duration)
        !wtime = omp_get_wtime()
        indices = argsort(eigvals(1, :))
        eigvectors_temp = eig_vector
        do count1 = 1, time_num
            do count2 = 1, sizeH
                energies(count1, count2) = eigvals(count1, indices(count2))
                eig_vector(count1, :, count2) = eigvectors_temp(count1, :, indices(count2))
            end do
        end do
        !wtime = omp_get_wtime () - wtime
        !write(6, *) 'Elapsed time:', sngl(wtime)

        ! Calculate inverse eigenvectors
        eig_vector_inv = inverse_real(eig_vector)

        ! Calculate Liouville space propagator with relaxation
        call liouville_propagator(time_num, time_step, electron_frequency, nuclear_frequency, microwave_amplitude, &
                t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, sizeH, sizeL, eig_vector, eig_vector_inv, energies, &
                propagator)

        ! Propagate density matrix stroboscopically, calculating polarisations
        call calculate_polarisation_rotor(time_num, time_num_prop, density_mat, propagator, sizeH, sizeL, &
                pol_i_z, pol_s1_z, pol_s2_z)

        ! Propagate density matrix for single rotor period, calculating polarisations
        call calculate_polarisation_sub_rotor(time_num, density_mat, propagator, sizeH, sizeL, &
                pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot)

    end subroutine main

    subroutine calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, temperature, hyperfine_coupling, &
                hyperfine_angles, orientation_ce_1, orientation_ce_2, electron_frequency, microwave_frequency, &
                nuclear_frequency, sizeH, sizeL, hamiltonian, density_mat)

        ! Calculate intrinsic Hilbert space Hamiltonian and density matrix from an ideal Hamiltonian
        ! Independent processes so number of OMP threads can take any value

        use omp_lib
        use interactions
        use iso_fortran_env
        use functions
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), parameter :: pi = 4._wp * atan(1._wp), Planck = 6.62607004E-34_wp, Boltzmann = 1.38064852E-23_wp

        integer, intent(in):: time_num, sizeH, sizeL
        real(wp), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_ce_1, orientation_ce_2
        real(wp), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(wp), intent(in) :: time_step, freq_rotor, temperature
        real(wp), dimension(time_num, sizeH, sizeH), intent(out) :: hamiltonian
        real(wp), intent(out) :: density_mat(sizeH, sizeH)

        integer :: count
        real(wp) :: c0_1, c1_1, c2_1, c3_1, c4_1, c0_2, c1_2, c2_2, c3_2, c4_2
        real(wp) :: hyperfine_zx, hyperfine_zz, ganisotropy_1, ganisotropy_2
        real(wp), dimension(2, 2) :: spin1_x, spin1_z, identity_size2

        real(wp), dimension(sizeH, sizeH) :: spin3_s1_z, spin3_s2_z, hyperfine_total, dipolar
        real(wp), dimension(sizeH, sizeH) :: spin3_i_x, spin3_i_z
        real(wp), dimension(sizeH, sizeH) :: hamiltonian_ideal, boltzmann_factors_mat
        real(wp), dimension(sizeH) :: boltzmann_factors

        real(wp) :: test4(4, 4), test8(8, 8), test8_2(sizeH, sizeH)

        ! Identity matrix
        identity_size2 = eye(2)

        ! Pauli matrices
        spin1_x = 0.5_wp * (reshape([0._wp, 1._wp, 1._wp, 0._wp], shape(spin1_x), order = [2, 1]))
        spin1_z = 0.5_wp * (reshape([ 1._wp, 0._wp, 0._wp, -1._wp], shape(spin1_z), order = [2, 1]))

        ! 4x4 spin matrices constructed using Kronecker products
        spin3_s1_z = kron_real(kron_real(spin1_z, identity_size2), identity_size2)
        spin3_s2_z = kron_real(kron_real(identity_size2, spin1_z), identity_size2)
        spin3_i_x = kron_real(kron_real(identity_size2, identity_size2), spin1_x)
        spin3_i_z = kron_real(kron_real(identity_size2, identity_size2), spin1_z)

        ! Calculate time independent electron g-anisotropy coefficients
        call anisotropy_coefficients(electron_frequency, gtensor, orientation_ce_1, c0_1, c1_1, c2_1, c3_1, c4_1)
        call anisotropy_coefficients(electron_frequency, gtensor, orientation_ce_2, c0_2, c1_2, c2_2, c3_2, c4_2)

        !$omp parallel do default(private) &
        !$omp& shared(c0_1, c1_1, c2_1, c3_1, c4_1, c0_2, c1_2, c2_2, c3_2, c4_2) &
        !$omp& shared(freq_rotor, time_step, hyperfine_angles, hyperfine_coupling) &
        !$omp& shared(hamiltonian, microwave_frequency, nuclear_frequency) &
        !$omp& shared(spin3_i_x, spin3_i_z, spin3_s1_z, spin3_s2_z)
        do count = 1, time_num

            ! Calculate time dependent hyperfine
            call hyperfine(hyperfine_coupling, hyperfine_angles, freq_rotor, (count - 1) * time_step, &
                    hyperfine_zz, hyperfine_zx)
            hyperfine_total = 2._wp * hyperfine_zx * matmul(spin3_i_x, spin3_s1_z) + &
                              2._wp * hyperfine_zz * matmul(spin3_i_z, spin3_s1_z)

            ! Calculate time dependent electron g-anisotropy
            call anisotropy(c0_1, c1_1, c2_1, c3_1, c4_1, freq_rotor, electron_frequency, (count - 1) * time_step, &
                    gtensor, ganisotropy_1)
            call anisotropy(c0_2, c1_2, c2_2, c3_2, c4_2, freq_rotor, electron_frequency, (count - 1) * time_step, &
                    gtensor, ganisotropy_2)

            ! Calculate time dependent dipolar between electron 1 and 2
            dipolar = 0

            ! Calculate time dependent hamiltonian
            hamiltonian(count, :, :) = (ganisotropy_1 - microwave_frequency) * spin3_s1_z + &
                                       (ganisotropy_2 - microwave_frequency) * spin3_s2_z + &
                                        nuclear_frequency * spin3_i_z + &
                                        hyperfine_total
        end do
        !$omp end parallel do

        ! Calculate idealised Hamiltonian (must calculate density matrix outside of rotating frame)
        hamiltonian_ideal = electron_frequency * (spin3_s1_z + spin3_s2_z) + nuclear_frequency * spin3_i_z

        ! Calculate initial Zeeman basis density matrix from Boltzmann factors
        do count = 1, sizeH
            boltzmann_factors(count) = exp(-(Planck * hamiltonian_ideal(count, count)) / (Boltzmann * temperature))
            boltzmann_factors_mat(count, count) = boltzmann_factors(count)
        end do
        density_mat = (1._wp / sum(boltzmann_factors)) * boltzmann_factors_mat

    end subroutine calculate_hamiltonian

    subroutine liouville_propagator(time_num, time_step, electron_frequency, nuclear_frequency, microwave_amplitude, &
                t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, sizeH, sizeL, eig_vector, eig_vector_inv, energies, &
                propagator)

        ! Calculate Louville space propagator with relaxation
        ! Independent processes so number of OMP threads can take any value

        use omp_lib
        use iso_fortran_env
        use functions
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), parameter :: i = (0, 1._wp)
        real(wp), parameter :: Planck = 6.62607004E-34_wp, Boltzmann = 1.38064852E-23_wp

        integer, intent(in):: time_num, sizeH, sizeL
        real(wp), dimension(:, :, :), intent(in) :: eig_vector, eig_vector_inv
        real(wp), dimension(:, :), intent(in) :: energies
        real(wp), intent(in) :: electron_frequency, nuclear_frequency, microwave_amplitude, time_step
        real(wp), intent(in) :: t1_nuc, t1_elec, t2_nuc, t2_elec, temperature
        complex(wp), intent(out) :: propagator(time_num, sizeL, sizeL)

        integer :: count, count2
        real(wp), dimension(2, 2) :: spin1_x, spin1_z, identity_size2
        complex(wp), dimension(2, 2) :: identity_size2_complex, spin1_y
        real(wp) :: p_e, p_n, gnp, gnm, gep, gem, boltzmann_elec, boltzmann_nuc

        real(wp), dimension(sizeL, sizeL) :: identity_size64, hamiltonian_liouville, relax_mat
        real(wp), dimension(sizeH, sizeH) :: spin3_s1_x, spin3_s2_x, spin3_s1_z, spin3_s2_z, total_hamiltonian
        real(wp), dimension(sizeH, sizeH) :: spin3_s1_p, spin3_s2_p, spin3_s1_m, spin3_s2_m, spin3_i_p, spin3_i_m
        real(wp), dimension(sizeH, sizeH) :: microwave_hamiltonian_init, microwave_hamiltonian, energy_mat
        real(wp), dimension(sizeH, sizeH) :: spin3_i_x, spin3_i_z, identity_size8

        complex(wp), dimension(sizeL, sizeL) :: eigvectors_liouville, eigvectors_inv_liouville
        complex(wp), dimension(sizeL, sizeL) :: liouvillian, mat_exp
        complex(wp), dimension(sizeH, sizeH) :: spin3_s1_y, spin3_s2_y, spin3_i_y

        ! Identity matrix
        identity_size2 = eye(2)
        identity_size2_complex = eye(2)
        identity_size8 = eye(8)
        identity_size64 = eye(64)

        ! Pauli matrices
        spin1_x = 0.5 * (reshape([0._wp, 1._wp, 1._wp, 0._wp], shape(spin1_x), order = [2, 1]))
        spin1_y = 0.5 * i * (reshape([0._wp, -1._wp, 1._wp, 0._wp], shape(spin1_y), order = [2, 1]))
        spin1_z = 0.5 * (reshape([1._wp, 0._wp, 0._wp, -1._wp], shape(spin1_z), order = [2, 1]))

        ! 4x4 matrices for S1 operator
        spin3_s1_x = kron_real(kron_real(spin1_x, identity_size2), identity_size2)
        spin3_s1_y = kron_complex(kron_complex(spin1_y, identity_size2_complex), identity_size2_complex)
        spin3_s1_z = kron_real(kron_real(spin1_z, identity_size2), identity_size2)
        spin3_s1_p = spin3_s1_x + real(i * spin3_s1_y)
        spin3_s1_m = spin3_s1_x - real(i * spin3_s1_y)

        ! 4x4 matrices for S2 operator
        spin3_s2_x = kron_real(kron_real(identity_size2, spin1_x), identity_size2)
        spin3_s2_y = kron_complex(kron_complex(identity_size2_complex, spin1_y), identity_size2_complex)
        spin3_s2_z = kron_real(kron_real(identity_size2, spin1_z), identity_size2)
        spin3_s2_p = spin3_s2_x + real(i * spin3_s2_y)
        spin3_s2_m = spin3_s2_x - real(i * spin3_s2_y)

        ! 4x4 matrices for I operator
        spin3_i_x = kron_real(kron_real(identity_size2, identity_size2), spin1_x)
        spin3_i_y = kron_complex(kron_complex(identity_size2_complex, identity_size2_complex), spin1_y)
        spin3_i_z = kron_real(kron_real(identity_size2, identity_size2), spin1_z)
        spin3_i_p = spin3_i_x + real(i * spin3_i_y)
        spin3_i_m = spin3_i_x - real(i * spin3_i_y)

        ! Calculate variables for origonal Liouville space relaxation
        p_e = tanh(0.5_wp * electron_frequency * (Planck / (Boltzmann * temperature)))
        p_n = tanh(0.5_wp * nuclear_frequency * (Planck / (Boltzmann * temperature)))
        gnp = 0.5_wp * (1._wp - p_n) * (1._wp / (1._wp * t1_nuc))
        gnm = 0.5_wp * (1._wp + p_n) * (1._wp / (1._wp * t1_nuc))
        gep = 0.5_wp * (1._wp - p_e) * (1._wp / (1._wp * t1_elec))
        gem = 0.5_wp * (1._wp + p_e) * (1._wp / (1._wp * t1_elec))

        ! Calculate variables for Mance Liouville space relaxation
        boltzmann_elec = electron_frequency * (Planck / (Boltzmann * temperature))
        boltzmann_nuc = nuclear_frequency * (Planck / (Boltzmann * temperature))

        ! Calculate initial microwave Hamiltonian
        microwave_hamiltonian_init = microwave_amplitude * (spin3_s1_x + spin3_s2_x)

        !$omp parallel do default(private) &
        !$omp& shared(eig_vector, eig_vector_inv, identity_size8, identity_size64, sizeL, sizeH) &
        !$omp& shared(t1_nuc, t1_elec, t2_elec, t2_nuc, boltzmann_elec, boltzmann_nuc) &
        !$omp& shared(spin3_s1_z, spin3_s1_x, spin3_s1_p, spin3_s1_m, spin3_s2_z, spin3_s2_p, spin3_s2_m, spin3_i_z) &
        !$omp& shared(spin3_s2_x, spin3_i_p, spin3_i_m, spin3_i_x) &
        !$omp& shared(gnp, gnm, gep, gem, microwave_hamiltonian_init, energies, propagator, time_step)
        do count = 1, time_num

            ! Transform microwave Hamiltonian into time dependent basis
            microwave_hamiltonian = matmul(eig_vector_inv(count, :, :), matmul(microwave_hamiltonian_init, &
                                            eig_vector(count, :, :)))

            ! Calculate total Hamiltonian
            energy_mat = 0._wp
            do count2 = 1, sizeH
                energy_mat(count2, count2) = energies(count, count2)
            end do

            ! Calculate total Hamiltonian
            total_hamiltonian = energy_mat + microwave_hamiltonian

            ! Transform Hilbert space Hamiltonian into Liouville space
            hamiltonian_liouville = kron_real(total_hamiltonian, identity_size8) - &
                    kron_real(identity_size8, transpose(total_hamiltonian))

            ! Calculate Louville space relaxation matrix using origonal theory
!            call calculate_relaxation_mat(eig_vector(count, :, :), eig_vector_inv(count, :, :), identity_size8, &
!                    identity_size64, sizeL, sizeH, spin3_s1_z, spin3_s1_p, spin3_s1_m, spin3_s2_z, spin3_s2_p, &
!                    spin3_s2_m, spin3_i_z, spin3_i_p, spin3_i_m, t2_elec, t2_nuc, gnp, gnm, gep, gem, relax_mat)

            ! Calculate Louville space relaxation matrix using Mance theory
            call calculate_relaxation_mat_mance(eig_vector(count, :, :), eig_vector_inv(count, :, :), sizeL, sizeH, &
                    spin3_s1_x, spin3_s1_z, spin3_s2_x, spin3_s2_z, spin3_i_x, spin3_i_z, &
                    t1_elec, t1_nuc, t2_elec, t2_nuc, boltzmann_elec, boltzmann_nuc, relax_mat)

            ! Calculate Liouville space eigenvectors
            eigvectors_liouville = kron_real(eig_vector(count, :, :), eig_vector(count, :, :))
            eigvectors_inv_liouville = kron_real(eig_vector_inv(count, :, :), eig_vector_inv(count, :, :))

            ! Calculate Liouville space propagator
            liouvillian = hamiltonian_liouville + i * relax_mat
            mat_exp = expm_complex(-i * liouvillian * time_step)
            propagator(count, :, :) = matmul(eigvectors_inv_liouville, matmul(mat_exp, eigvectors_liouville))

        end do
        !$omp end parallel do

    end subroutine liouville_propagator

    subroutine calculate_relaxation_mat(eig_vector, eig_vector_inv, identity_size8, identity_size64, sizeL, sizeH, &
            spin3_s1_z, spin3_s1_p, spin3_s1_m, spin3_s2_z, spin3_s2_p, spin3_s2_m, spin3_i_z, spin3_i_p, spin3_i_m, &
            t2_elec, t2_nuc, gnp, gnm, gep, gem, relax_mat)

        ! Calculate Louville space relaxation matrix

        use functions
        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        integer, intent(in) :: sizeH, sizeL
        real(wp), dimension(:, :), intent(in) :: eig_vector, eig_vector_inv
        real(wp), dimension(:, :), intent(in) :: spin3_s1_z, spin3_s1_p, spin3_s1_m, spin3_s2_z, spin3_s2_p, spin3_s2_m
        real(wp), dimension(:, :), intent(in) :: spin3_i_z, spin3_i_p, spin3_i_m, identity_size8
        real(wp), dimension(:, :), intent(in) :: identity_size64
        real(wp), intent(in) :: t2_nuc, t2_elec
        real(wp), intent(in) :: gnp, gnm, gep, gem

        real(wp), dimension(sizeL, sizeL), intent(out) :: relax_mat

        real(wp), dimension(sizeL, sizeL) :: spin3_i_p_tl, spin3_i_m_tl, spin3_s1_p_tl, spin3_s1_m_tl, spin3_s2_p_tl
        real(wp), dimension(sizeL, sizeL) :: spin3_s2_m_tl, relax_t2_s1, relax_t2_s2, relax_t2_i1, relax_t1
        real(wp), dimension(sizeH, sizeH) :: spin3_s1_z_t, spin3_s1_p_t, spin3_s1_m_t, spin3_i_z_t, spin3_i_p_t
        real(wp), dimension(sizeH, sizeH) :: spin3_i_m_t, spin3_s2_z_t, spin3_s2_p_t, spin3_s2_m_t

        ! Transform spin matrices into time dependent Hilbert space basis
        spin3_s1_z_t = matmul(eig_vector_inv, matmul(spin3_s1_z, eig_vector))
        spin3_s1_p_t = matmul(eig_vector_inv, matmul(spin3_s1_p, eig_vector))
        spin3_s1_m_t = matmul(eig_vector_inv, matmul(spin3_s1_m, eig_vector))

        spin3_s2_z_t = matmul(eig_vector_inv, matmul(spin3_s2_z, eig_vector))
        spin3_s2_p_t = matmul(eig_vector_inv, matmul(spin3_s2_p, eig_vector))
        spin3_s2_m_t = matmul(eig_vector_inv, matmul(spin3_s2_m, eig_vector))

        spin3_i_z_t = matmul(eig_vector_inv, matmul(spin3_i_z, eig_vector))
        spin3_i_p_t = matmul(eig_vector_inv, matmul(spin3_i_p, eig_vector))
        spin3_i_m_t = matmul(eig_vector_inv, matmul(spin3_i_m, eig_vector))

        ! Transform spin matrices into time dependent Liouville space basis
        spin3_i_p_tl = kron_real(spin3_i_p_t, transpose(spin3_i_m_t)) - 0.5_wp * identity_size64 + 0.5_wp * (&
                    kron_real(spin3_i_z_t, identity_size8) + kron_real(identity_size8, transpose(spin3_i_z_t)))

        spin3_i_m_tl = kron_real(spin3_i_m_t, transpose(spin3_i_p_t)) - 0.5_wp * identity_size64 - 0.5_wp * (&
                kron_real(spin3_i_z_t, identity_size8) + kron_real(identity_size8, transpose(spin3_i_z_t)))

        spin3_s1_p_tl = kron_real(spin3_s1_p_t, transpose(spin3_s1_m_t)) - 0.5_wp * identity_size64 + 0.5_wp * (&
                kron_real(spin3_s1_z_t, identity_size8) + kron_real(identity_size8, transpose(spin3_s1_z_t)))

        spin3_s1_m_tl = kron_real(spin3_s1_m_t, transpose(spin3_s1_p_t)) - 0.5_wp * identity_size64 - 0.5_wp * (&
                kron_real(spin3_s1_z_t, identity_size8) + kron_real(identity_size8, transpose(spin3_s1_z_t)))

        spin3_s2_p_tl = kron_real(spin3_s2_p_t, transpose(spin3_s2_m_t)) - 0.5_wp * identity_size64 + 0.5_wp * (&
                kron_real(spin3_s2_z_t, identity_size8) + kron_real(identity_size8, transpose(spin3_s2_z_t)))

        spin3_s2_m_tl = kron_real(spin3_s2_m_t, transpose(spin3_s2_p_t)) - 0.5_wp * identity_size64 - 0.5_wp * (&
                kron_real(spin3_s2_z_t, identity_size8) + kron_real(identity_size8, transpose(spin3_s2_z_t)))

        ! Calculate time dependent Liouville space relaxation matrix
        relax_t2_s1 = (2._wp / t2_elec) * (kron_real(spin3_s1_z_t, transpose(spin3_s1_z_t)) - &
                        0.5_wp * 0.5_wp * identity_size64)
        relax_t2_s2 = (2._wp / t2_elec) * (kron_real(spin3_s2_z_t, transpose(spin3_s2_z_t)) - &
                        0.5_wp * 0.5_wp * identity_size64)
        relax_t2_i1 = (2._wp / t2_nuc) * (kron_real(spin3_i_z_t, transpose(spin3_i_z_t)) - &
                        0.5_wp * 0.5_wp * identity_size64)
        relax_t1 = gep * spin3_s1_p_tl + gem * spin3_s1_m_tl + &
                    gep * spin3_s2_p_tl + gem * spin3_s2_m_tl + &
                    gnp * spin3_i_p_tl + gnm * spin3_i_m_tl
        relax_mat = relax_t2_s1 + relax_t2_s2 + relax_t2_i1 + relax_t1

    end subroutine calculate_relaxation_mat

    subroutine calculate_relaxation_mat_mance(eig_vector, eig_vector_inv, sizeL, sizeH, &
            spin3_s1_x, spin3_s1_z, spin3_s2_x, spin3_s2_z, spin3_i_x, spin3_i_z, &
            t1_elec, t1_nuc, t2_elec, t2_nuc, boltzmann_elec, boltzmann_nuc, relax_mat)

        ! Calculate Louville space relaxation matrix using Mance theory

        use functions
        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        integer, intent(in) :: sizeH, sizeL
        real(wp), dimension(:, :), intent(in) :: eig_vector, eig_vector_inv, spin3_i_z
        real(wp), dimension(:, :), intent(in) :: spin3_s1_x, spin3_s1_z, spin3_s2_x, spin3_s2_z, spin3_i_x
        real(wp), intent(in) :: t1_elec, t1_nuc, t2_nuc, t2_elec, boltzmann_elec, boltzmann_nuc

        real(wp), dimension(sizeL, sizeL), intent(out) :: relax_mat

        real(wp), dimension(sizeH, sizeH) :: spin3_s1_x_t, spin3_s1_z_t, spin3_s2_x_t, spin3_s2_z_t
        real(wp), dimension(sizeH, sizeH) :: spin3_i_x_t, spin3_i_z_t, relax_values_t1
        real(wp), dimension(sizeL, sizeL) :: relax_mat_t2, relax_mat_t1
        real(wp) :: relax_values_t2(sizeH), diff_elec, diff_nuc, boltzmann_factor
        integer :: count, count2, count3, count4

        ! Transform spin matrices into time dependent Hilbert space basis
        spin3_s1_x_t = matmul(eig_vector_inv, matmul(spin3_s1_x, eig_vector))
        spin3_s1_z_t = matmul(eig_vector_inv, matmul(spin3_s1_z, eig_vector))
        spin3_s2_x_t = matmul(eig_vector_inv, matmul(spin3_s2_x, eig_vector))
        spin3_s2_z_t = matmul(eig_vector_inv, matmul(spin3_s2_z, eig_vector))
        spin3_i_x_t = matmul(eig_vector_inv, matmul(spin3_i_x, eig_vector))
        spin3_i_z_t = matmul(eig_vector_inv, matmul(spin3_i_z, eig_vector))

        ! Calculate T1 relaxation matrix using Boltzmann factors
        do count = 1, sizeH
            do count2 = 1, sizeH

                ! Calculate Boltzmann factors
                diff_elec = spin3_s1_z_t(count, count) - spin3_s1_z_t(count2, count2) + &
                            spin3_s2_z_t(count, count) - spin3_s2_z_t(count2, count2)
                diff_nuc = spin3_i_z_t(count, count) - spin3_i_z_t(count2, count2)
                boltzmann_factor = exp(-diff_elec * boltzmann_elec / 2._wp - diff_nuc * boltzmann_nuc / 2._wp) / &
                                   (exp(diff_elec * boltzmann_elec / 2._wp + diff_nuc * boltzmann_nuc / 2._wp) + &
                                    exp(-diff_elec * boltzmann_elec / 2._wp - diff_nuc * boltzmann_nuc / 2._wp))

                ! Calculate initial relaxation matrix
                relax_values_t1(count, count2) = &
                        (1._wp / t1_elec) * (spin3_s1_x_t(count, count2) * spin3_s1_x_t(count2, count) + &
                                             spin3_s2_x_t(count, count2) * spin3_s2_x_t(count2, count) + &
                                             spin3_s1_z_t(count, count2) * spin3_s1_z_t(count2, count) + &
                                             spin3_s2_z_t(count, count2) * spin3_s2_z_t(count2, count)) + &
                        (1._wp / t1_nuc) * (spin3_i_x_t(count, count2) * spin3_i_x_t(count2, count) + &
                                            spin3_i_z_t(count, count2) * spin3_i_z_t(count2, count))

                relax_values_t1(count, count2) = relax_values_t1(count, count2) * boltzmann_factor
            end do
        end do

        ! Set diagonals equal to zero
        do count = 1, sizeH
            relax_values_t1(count, count) = 0._wp
        end do

        ! Set diagonal elements
        do count = 1, sizeH
            do count2 = 1, sizeH
                if (abs(count - count2) > 0) then
                    relax_values_t1(count, count) = relax_values_t1(count, count) - relax_values_t1(count2, count)
                end if
            end do
        end do

        ! Transform to 16 by 16 matrix
        relax_mat_t1 = 0
        do count = 1, sizeH
            do count2 = 1, sizeH
                count3 = (count-1) * sizeH + count
                count4 = (count2-1) * sizeH + count2
                relax_mat_t1(count3, count4) = relax_values_t1(count, count2)
            end do
        end do

        ! Calculate T2 relaxation matrix using circular shift of spin matrices
        relax_values_t2 = 0
        do count = 2, sizeH
            relax_values_t2(count) = count-1
            relax_values_t2(count) = &
                ((abs(spin3_s1_z_t(count, count) - spin3_s1_z_t(count - 1, count - 1))) ** 2._wp + &
                 (abs(spin3_s2_z_t(count, count) - spin3_s2_z_t(count - 1, count - 1))) ** 2._wp) * &
                        (1._wp / t2_elec) + &
                ((abs(spin3_i_z_t(count, count) - spin3_i_z_t(count - 1, count - 1))) ** 2._wp) * (1._wp / t2_nuc)
        end do

!        write(6, *) 'relax_values_t2 before', relax_values_t2
        !relax_values_t2 = cshift(relax_values_t2, shift = -1)
        !write(6, *) 'relax_values_t2 after', relax_values_t2

        ! Transform to 16 by 16 matrix
        relax_mat_t2 = 0
        do count = 1, sizeH
            do count2 = 1, sizeH
                relax_mat_t2(count2 + sizeH * (count - 1), count2 + sizeH * (count - 1)) = relax_values_t2(count2)
            end do
            relax_values_t2 = cshift(relax_values_t2, shift = -1)
!            write(6, *) 'count', count
!            write(6, *) 'relax_values_t2 after', relax_values_t2
        end do

        ! Relaxation matrix as sum of t1 and t2 matrices
        relax_mat = -1._wp * relax_mat_t2 + relax_mat_t1
        !relax_mat = relax_mat_t1

!        write(6, *) 'relax_mat', relax_mat(9:12, 9:12)

!        write(6, *) 'relax_mat_t2(1, 1)', relax_mat_t2(1, 1)
!        write(6, *) 'relax_mat_t2(1, 1)', relax_mat_t2(1, 2)
!        write(6, *) 'relax_mat_t2(2, 2)', relax_mat_t2(2, 2)
!        write(6, *) 'relax_mat_t2(9, 9)', relax_mat_t2(9, 9)
!        write(6, *) 'relax_mat_t2(10, 10)', relax_mat_t2(10, 10)
!        write(6, *) 'relax_mat_t2(64, 64)', relax_mat_t2(63, 63)
!        write(6, *) 'relax_mat_t2(64, 64)', relax_mat_t2(64, 64)


!        write(6, *) 'relax_mat_t1(1, 1)', relax_mat_t1(1, 1)
!        write(6, *) 'relax_mat_t1(2, 2)', relax_mat_t1(2, 2)

    end subroutine calculate_relaxation_mat_mance

    subroutine calculate_polarisation_rotor(time_num, time_num_prop, density_mat, propagator, sizeH, sizeL, &
                pol_i_z, pol_s1_z, pol_s2_z)

        ! Calculate polarisation stroboscopically across multiple rotor periods
        ! Iterative process so number of OMP threads must equal 1

        use omp_lib
        use functions
        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        integer, intent(in) :: time_num_prop, time_num, sizeH, sizeL
        complex(wp), dimension(:, :, :), intent(in) :: propagator
        real(wp), dimension(:, :), intent(in) :: density_mat
        real(wp), dimension(time_num_prop), intent(out) :: pol_i_z, pol_s1_z, pol_s2_z

        integer :: count
        complex(wp), dimension(2, 2) :: spin1_z, identity_size2

        complex(wp), dimension(sizeL, 1) :: density_mat_liouville, temp
        complex(wp), dimension(sizeL, sizeL) :: identity_size64, propagator_strobe
        complex(wp), dimension(sizeH, sizeH) :: spin3_s1_z, spin3_s2_z, spin3_i_z, density_mat_time

        ! Calculate matrices specific to polarisation calculation
        identity_size2 = eye(2)
        identity_size64 = eye(64)

        ! Calculate matrices specific to polarisation calculation
        spin1_z = 0.5_wp * (reshape([1._wp, 0._wp,  0._wp, -1._wp], shape(spin1_z), order = [2, 1]))
        spin3_s1_z = kron_complex(kron_complex(spin1_z, identity_size2), identity_size2)
        spin3_s2_z = kron_complex(kron_complex(identity_size2, spin1_z), identity_size2)
        spin3_i_z = kron_complex(kron_complex(identity_size2, identity_size2), spin1_z)

        density_mat_time = density_mat
        propagator_strobe = identity_size64

        ! Calculate stroboscopic propagator (product of all operators within rotor period)
        do count = 1, time_num
            propagator_strobe = matmul(propagator_strobe, propagator(count, :, :))

        end do

        ! Propogate density matrix using stroboscopic propagator
        do count = 1, time_num_prop

            ! Calculate electronic and nuclear polarisation
            pol_i_z(count) = real(trace_complex(matmul(density_mat_time, spin3_i_z)))
            pol_s1_z(count) = real(trace_complex(matmul(density_mat_time, spin3_s1_z)))
            pol_s2_z(count) = real(trace_complex(matmul(density_mat_time, spin3_s2_z)))

            ! Transform density matrix (2^N x 2^N to 4^N x 1)
            density_mat_liouville = reshape(density_mat_time, [sizeL, 1])

            ! Propagate density matrix
            temp = matmul(propagator_strobe, density_mat_liouville)

            ! Transform density matrix (4^N x 1 to 2^N x 2^N)
            density_mat_time = reshape(temp, [sizeH, sizeH])

        end do

    end subroutine calculate_polarisation_rotor

    subroutine calculate_polarisation_sub_rotor(time_num, density_mat, propagator, sizeH, sizeL, &
                pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot)

        ! Calculate sub rotor polarisation
        ! Iterative process so number of OMP threads must equal 1

        use omp_lib
        use functions
        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        integer, intent(in) :: time_num, sizeH, sizeL
        complex(wp), dimension(:, :, :), intent(in) :: propagator
        real(wp), dimension(:, :), intent(in) :: density_mat
        real(wp), dimension(time_num), intent(out) :: pol_i_z_rot, pol_s1_z_rot, pol_s2_z_rot

        integer :: count
        complex(wp), dimension(2, 2) :: spin1_z, identity_size2

        complex(wp), dimension(sizeL, 1) :: density_mat_liouville, temp
        complex(wp), dimension(sizeH, sizeH) :: spin3_s1_z, spin3_s2_z, spin3_i_z, density_mat_time

        ! Calculate matrices specific to polarisation calculation
        identity_size2 = transpose(reshape([1._wp, 0._wp, 0._wp, 1._wp], shape(identity_size2)))
        spin1_z = 0.5_wp * (reshape([1._wp, 0._wp, 0._wp, -1._wp], shape(spin1_z), order = [2, 1]))
        spin3_s1_z = kron_complex(kron_complex(spin1_z, identity_size2), identity_size2)
        spin3_s2_z = kron_complex(kron_complex(identity_size2, spin1_z), identity_size2)
        spin3_i_z = kron_complex(kron_complex(identity_size2, identity_size2), spin1_z)

        density_mat_time = density_mat

        ! Propogate density matrix using stroboscopic propagator
        do count = 1, time_num

            ! Calculate electronic and nuclear polarisation
            pol_i_z_rot(count) = real(trace_complex(matmul(density_mat_time, spin3_i_z)))
            pol_s1_z_rot(count) = real(trace_complex(matmul(density_mat_time, spin3_s1_z)))
            pol_s2_z_rot(count) = real(trace_complex(matmul(density_mat_time, spin3_s2_z)))

            ! Transform density matrix (2^N x 2^N to 4^N x 1)
            density_mat_liouville = reshape(density_mat_time, [sizeL, 1])

            ! Propagate density matrix
            temp = matmul(propagator(count, :, :), density_mat_liouville)

            ! Transform density matrix (4^N x 1 to 2^N x 2^N)
            density_mat_time = reshape(temp, [sizeH, sizeH])

        end do

    end subroutine calculate_polarisation_sub_rotor

end module
