module solid_effect_dynamics

    ! This module calculates dynamics for solid effect MAS DNP NMR.

contains

    subroutine main(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, &
                orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, microwave_amplitude, &
                t1_nuc,  t1_elec, t2_nuc, t2_elec, temperature, time_num_prop, &
                pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies)

        ! Call required functions to calculate dynamics for solid effect MAS DNP NMR

         use omp_lib
         use functions
        implicit none

        integer, parameter :: sizeH = 2 ** (2), sizeL = 4 ** (2) ! Change (value) to change number of spins

        integer, intent(in):: time_num, time_num_prop
        real(kind=8), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_se
        real(kind=8), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(kind=8), intent(in) :: time_step, freq_rotor, microwave_amplitude
        real(kind=8), intent(in) :: t1_nuc, t1_elec, t2_nuc, t2_elec, temperature

        real(kind=8), dimension(time_num_prop), intent(out) :: pol_i_z, pol_s_z
        real(kind=8), dimension(time_num), intent(out) :: pol_i_z_rot, pol_s_z_rot
        real(kind=8), dimension(time_num, sizeH), intent(out) :: energies

        complex(kind=8), dimension(time_num, sizeL, sizeL)  :: propagator
        real(kind=8), dimension(sizeH, sizeH) :: density_mat
        real(kind=8) :: eig_vector(time_num, sizeH,sizeH), eig_vector_inv(time_num, sizeH,sizeH)
        real(kind=8), dimension(time_num, sizeH, sizeH) :: hamiltonian

        ! Construct intrinsic Hilbert space Hamiltonian
        call calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, temperature, hyperfine_coupling, &
                hyperfine_angles, orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, &
                sizeH, sizeL, hamiltonian, density_mat)

        ! Calculate eigenvalues and eigenvectors of intrinsic Hamiltonian
        call eig_real(hamiltonian, energies, eig_vector)
        eig_vector_inv = inverse_real(eig_vector)

        ! Calculate Liouville space propagator with relaxation
        call liouville_propagator(time_num, time_step, electron_frequency, nuclear_frequency, microwave_amplitude, &
                t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, sizeH, sizeL, eig_vector, eig_vector_inv, energies, &
                propagator)

        ! Propagate density matrix stroboscopically, calculating polarisations
        call calculate_polarisation_rotor(time_num, time_num_prop, density_mat, propagator, sizeH, sizeL, &
                pol_i_z, pol_s_z)

        ! Propagate density matrix for single rotor period, calculating polarisations
        call calculate_polarisation_sub_rotor(time_num, density_mat, propagator, sizeH, sizeL, &
                pol_i_z_rot, pol_s_z_rot)

    end subroutine main

    subroutine calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, temperature, hyperfine_coupling, &
                hyperfine_angles, orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, &
                sizeH, sizeL, hamiltonian, density_mat)

        ! Calculate intrinsic Hilbert space Hamiltonian and density matrix from an ideal Hamiltonian
        ! Independent processes so number of OMP threads can take any value

        use omp_lib
        use interactions
        use functions
        implicit none

        real(kind=8), parameter :: PI = 4.D0 * DATAN(1.D0), Planck = 6.62607004D-34, Boltzmann = 1.38064852D-23

        integer, intent(in):: time_num, sizeH, sizeL
        real(kind=8), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_se
        real(kind=8), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(kind=8), intent(in) :: time_step, freq_rotor, temperature
        real(kind=8), dimension(time_num, sizeH, sizeH), intent(out) :: hamiltonian
        real(kind=8), intent(out) :: density_mat(sizeH, sizeH)

        integer :: count
        real(kind=8) :: c0, c1, c2, c3, c4
        real(kind=8) :: hyperfine_zx, hyperfine_zz, ganisotropy
        real(kind=8), dimension(2, 2) :: spin_x, spin_z, identity_spin1

        real(kind=8), dimension(sizeH, sizeH) :: spin2_s_x, spin2_s_z, hyperfine_total
        real(kind=8), dimension(sizeH, sizeH) :: spin2_i_x, spin2_i_z
        real(kind=8), dimension(sizeH, sizeH) :: hamiltonian_ideal, boltzmann_factors_mat
        real(kind=8), dimension(sizeH) :: boltzmann_factors

        ! Identity matrix
        identity_spin1 = transpose(reshape((/ 1.D0, 0.D0, 0.D0, 1.D0/), shape(identity_spin1)))

        ! Pauli matrices
        spin_x = 0.5D0 * (reshape((/ 0.D0, 1.D0, 1.D0, 0.D0/), shape(spin_x), order = (/2, 1/)))
        spin_z = 0.5D0 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron_real(spin_x, identity_spin1)
        spin2_s_z = kron_real(spin_z, identity_spin1)

        ! 4x4 matrices for I operator
        spin2_i_x = kron_real(identity_spin1, spin_x)
        spin2_i_z = kron_real(identity_spin1, spin_z)

        ! Calculate time independent electron g-anisotropy coefficients
        call anisotropy_coefficients(electron_frequency, gtensor, orientation_se, c0, c1, c2, c3, c4)

        !$omp parallel do default(private) &
        !$omp& shared(c0, c1, c2, c3, c4, freq_rotor, time_step, hyperfine_angles, hyperfine_coupling) &
        !$omp& shared(hamiltonian, microwave_frequency, nuclear_frequency) &
        !$omp& shared(spin2_i_x, spin2_i_z, spin2_s_z)
        do count = 1, time_num

            ! Calculate time dependent hyperfine
            call hyperfine(hyperfine_coupling, hyperfine_angles, freq_rotor, (count - 1) * time_step, &
                    hyperfine_zz, hyperfine_zx)
            hyperfine_total = 2.D0 * hyperfine_zx * MATMUL(spin2_i_x, spin2_s_z) + &
                              2.D0 * hyperfine_zz * MATMUL(spin2_i_z, spin2_s_z)

            ! Calculate time dependent electron g-anisotropy
            call anisotropy(c0, c1, c2, c3, c4, freq_rotor, electron_frequency, (count - 1) * time_step, &
                    gtensor, ganisotropy)

            ! Calculate time dependent hamiltonian
            hamiltonian(count, :, :) = (ganisotropy - microwave_frequency) * spin2_s_z + &
                                        nuclear_frequency * spin2_i_z + &
                                        hyperfine_total

        end do
        !$omp end parallel do

        ! Calculate idealised Hamiltonian (must calculate density matrix outside of rotating frame)
        hamiltonian_ideal = electron_frequency * spin2_s_z + nuclear_frequency * spin2_i_z

        ! Calculate initial Zeeman basis density matrix from Boltzmann factors
        !$omp parallel do default(private) &
        !$omp& shared(boltzmann_factors, boltzmann_factors_mat, hamiltonian_ideal, temperature)
        do count = 1, sizeH
            boltzmann_factors(count) = exp(-(Planck * hamiltonian_ideal(count, count)) / (Boltzmann * temperature))
            boltzmann_factors_mat(count, count) = boltzmann_factors(count)
        end do
        !$omp end parallel do
        density_mat = (1/sum(boltzmann_factors)) * boltzmann_factors_mat


    end subroutine calculate_hamiltonian

    subroutine liouville_propagator(time_num, time_step, electron_frequency, nuclear_frequency, microwave_amplitude, &
                t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, sizeH, sizeL, eig_vector, eig_vector_inv, energies, &
                propagator)

        ! Calculate Louville space propagator with relaxation
        ! Independent processes so number of OMP threads can take any value

        use omp_lib
        use functions
        implicit none

        complex(kind=8), parameter :: i = (0, 1.D0)
        real(kind=8), parameter :: Planck = 6.62607004D-34, Boltzmann = 1.38064852D-23

        integer, intent(in):: time_num, sizeH, sizeL
        real(kind=8), dimension(time_num, sizeH, sizeH), intent(in) :: eig_vector, eig_vector_inv
        real(kind=8), dimension(time_num, sizeH), intent(in) :: energies
        real(kind=8), intent(in) :: electron_frequency, nuclear_frequency, microwave_amplitude, time_step
        real(kind=8), intent(in) :: t1_nuc, t1_elec, t2_nuc, t2_elec, temperature
        complex(kind=8), intent(out) :: propagator(time_num, sizeL, sizeL)

        integer :: count, count2
        real(kind=8), dimension(2, 2) :: spin_x, spin_z, identity_size2
        complex(kind=8), dimension(2, 2) :: identity_size2_complex, spin_y
        real(kind=8) :: p_e, p_n, gnp, gnm, gep, gem

        real(kind=8), dimension(sizeL, sizeL) :: identity_size16, hamiltonian_liouville, relax_mat
        real(kind=8), dimension(sizeH, sizeH) :: spin2_s_x, spin2_s_z, total_hamiltonian
        real(kind=8), dimension(sizeH, sizeH) :: spin2_s_p, spin2_s_m, spin2_i_p, spin2_i_m
        real(kind=8), dimension(sizeH, sizeH) :: microwave_hamiltonian_init, microwave_hamiltonian, energy_mat
        real(kind=8), dimension(sizeH, sizeH) :: spin2_i_x, spin2_i_z, identity_size4

        complex(kind=8), dimension(sizeL, sizeL) :: eigvectors_liouville, eigvectors_inv_liouville
        complex(kind=8), dimension(sizeL, sizeL) :: liouvillian, mat_exp
        complex(kind=8), dimension(sizeH, sizeH) :: identity_size4_complex, spin2_s_y, spin2_i_y

        ! Identity matrix
        identity_size2 = transpose(reshape((/ 1.D0, 0.D0, 0.D0, 1.D0/), shape(identity_size2)))
        identity_size2_complex = transpose(reshape((/ 1.D0, 0.D0, 0.D0, 1.D0/), shape(identity_size2)))
        identity_size4 = kron_real(identity_size2, identity_size2)
        identity_size4_complex = kron_complex(identity_size2_complex, identity_size2_complex)
        identity_size16 = kron_real(kron_real(identity_size2, identity_size2), &
                kron_real(identity_size2, identity_size2))

        ! Pauli matrices
        spin_x = 0.5 * (reshape((/ 0.D0, 1.D0, 1.D0, 0.D0/), shape(spin_x), order = (/2, 1/)))
        spin_y = 0.5 * i * (reshape((/ 0.D0, -1.D0, 1.D0, 0.D0/), shape(spin_y), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron_real(spin_x, identity_size2)
        spin2_s_y = kron_complex(spin_y, identity_size2_complex)
        spin2_s_z = kron_real(spin_z, identity_size2)
        spin2_s_p = spin2_s_x + real(i * spin2_s_y)
        spin2_s_m = spin2_s_x - real(i * spin2_s_y)

        ! 4x4 matrices for I operator
        spin2_i_x = kron_real(identity_size2, spin_x)
        spin2_i_y = kron_complex(identity_size2_complex, spin_y)
        spin2_i_z = kron_real(identity_size2, spin_z)
        spin2_i_p = spin2_i_x + real(i * spin2_s_y)
        spin2_i_m = spin2_i_x - real(i * spin2_s_y)

        ! Calculate variables for Liouville space relaxation
        p_e = tanh(0.5D0 * electron_frequency * (Planck / (Boltzmann * temperature)))
        p_n = tanh(0.5D0 * nuclear_frequency * (Planck / (Boltzmann * temperature)))
        gnp = 0.5D0 * (1.D0 - p_n) * (1.D0 / (1.D0 * t1_nuc))
        gnm = 0.5D0 * (1.D0 + p_n) * (1.D0 / (1.D0 * t1_nuc))
        gep = 0.5D0 * (1.D0 - p_e) * (1.D0 / (1.D0 * t1_elec))
        gem = 0.5D0 * (1.D0 + p_e) * (1.D0 / (1.D0 * t1_elec))

        ! Calculate initial microwave Hamiltonian
        microwave_hamiltonian_init = microwave_amplitude * spin2_s_x
        
        !$omp parallel do default(private) &
        !$omp& shared(eig_vector, eig_vector_inv, identity_size4, identity_size16, sizeL, sizeH) &
        !$omp& shared(spin2_s_z, spin2_s_p, spin2_s_m, spin2_i_z, spin2_i_p, spin2_i_m, t2_elec, t2_nuc) &
        !$omp& shared(gnp, gnm, gep, gem, microwave_hamiltonian_init, energies, propagator, time_step)
        do count = 1, time_num

            ! Transform microwave Hamiltonian into time dependent basis
            microwave_hamiltonian = matmul(eig_vector_inv(count, :, :), matmul(microwave_hamiltonian_init, &
                    eig_vector(count, :, :)))

            ! Calculate total Hamiltonian
            energy_mat = identity_size4
            do count2 = 1, sizeH
                energy_mat(count2, count2) = energies(count, count2)
            end do

            ! Calculate total Hamiltonian
            total_hamiltonian = energy_mat + microwave_hamiltonian

            ! Transform Hilbert space Hamiltonian into Liouville space
            hamiltonian_liouville = kron_real(total_hamiltonian, identity_size4) - &
                    kron_real(identity_size4, transpose(total_hamiltonian))

            ! Calculate Louville space relaxation matrix
            call calculate_relaxation_mat(eig_vector(count, :, :), eig_vector_inv(count, :, :), identity_size4, &
                    identity_size16, sizeL, sizeH, spin2_s_z, spin2_s_p, spin2_s_m, spin2_i_z, spin2_i_p, spin2_i_m, &
                    t2_elec, t2_nuc, gnp, gnm, gep, gem, relax_mat)

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

    subroutine calculate_relaxation_mat(eig_vector, eig_vector_inv, identity_size4, identity_size16, sizeL, sizeH, &
            spin2_s_z, spin2_s_p, spin2_s_m, spin2_i_z, spin2_i_p, spin2_i_m, t2_elec, t2_nuc, &
            gnp, gnm, gep, gem, relax_mat)

        ! Calculate Louville space relaxation matrix

        use functions
        implicit none

        integer, intent(in) :: sizeH, sizeL
        real(kind = 8), dimension(sizeH, sizeH), intent(in) :: eig_vector, eig_vector_inv
        real(kind = 8), dimension(sizeH, sizeH), intent(in) :: spin2_s_z, spin2_i_z, identity_size4
        real(kind = 8), dimension(sizeH, sizeH), intent(in) :: spin2_s_p, spin2_s_m, spin2_i_p, spin2_i_m
        real(kind=8), dimension(sizeL, sizeL), intent(in) :: identity_size16
        real(kind = 8), intent(in) :: t2_nuc, t2_elec
        real(kind = 8), intent(in) :: gnp, gnm, gep, gem

        real(kind = 8), dimension(sizeL, sizeL), intent(out) :: relax_mat

        real(kind = 8), dimension(sizeL, sizeL) :: spin2_i_p_tl, spin2_i_m_tl, spin2_s_p_tl, spin2_s_m_tl
        real(kind = 8), dimension(sizeL, sizeL) :: relax_t2_elec, relax_t2_nuc, relax_t1
        real(kind = 8), dimension(sizeH, sizeH) :: spin2_s_z_t, spin2_s_p_t, spin2_s_m_t, spin2_i_z_t, spin2_i_p_t
        real(kind = 8), dimension(sizeH, sizeH) :: spin2_i_m_t

        ! Transform spin matrices into time dependent Hilbert space basis
        spin2_s_z_t = matmul(eig_vector_inv, matmul(spin2_s_z, eig_vector))
        spin2_s_p_t = matmul(eig_vector_inv, matmul(spin2_s_p, eig_vector))
        spin2_s_m_t = matmul(eig_vector_inv, matmul(spin2_s_m, eig_vector))
        spin2_i_z_t = matmul(eig_vector_inv, matmul(spin2_i_z, eig_vector))
        spin2_i_p_t = matmul(eig_vector_inv, matmul(spin2_i_p, eig_vector))
        spin2_i_m_t = matmul(eig_vector_inv, matmul(spin2_i_m, eig_vector))

        ! Transform spin matrices into time dependent Liouville space basis
        spin2_i_p_tl = kron_real(spin2_i_p_t, transpose(spin2_i_m_t)) - 0.5D0 * identity_size16 + 0.5D0 * (&
                kron_real(spin2_i_z_t, identity_size4) + &
                        kron_real(identity_size4, transpose(spin2_i_z_t)))

        spin2_i_m_tl = kron_real(spin2_i_m_t, transpose(spin2_i_p_t)) - 0.5D0 * identity_size16 - 0.5D0 * (&
                kron_real(spin2_i_z_t, identity_size4) + &
                        kron_real(identity_size4, transpose(spin2_i_z_t)))

        spin2_s_p_tl = kron_real(spin2_s_p_t, transpose(spin2_s_m_t)) - 0.5D0 * identity_size16 + 0.5D0 * (&
                kron_real(spin2_s_z_t, identity_size4) + &
                        kron_real(identity_size4, transpose(spin2_s_z_t)))

        spin2_s_m_tl = kron_real(spin2_s_m_t, transpose(spin2_s_p_t)) - 0.5D0 * identity_size16 - 0.5D0 * (&
                kron_real(spin2_s_z_t, identity_size4) + &
                        kron_real(identity_size4, transpose(spin2_s_z_t)))

        ! Calculate time dependent Liouville space relaxation matrix
        relax_t2_elec = (1.D0 / t2_elec) * (kron_real(spin2_s_z_t, transpose(spin2_s_z_t)) - &
                0.5D0 * 0.5D0 * identity_size16)
        relax_t2_nuc = (1.D0 / t2_nuc) * (kron_real(spin2_i_z_t, transpose(spin2_i_z_t)) - &
                0.5D0 * 0.5D0 * identity_size16)
        relax_t1 = gep * spin2_s_p_tl + gem * spin2_s_m_tl + gnp * spin2_i_p_tl + gnm * spin2_i_m_tl
        relax_mat = relax_t2_elec + relax_t2_nuc + relax_t1

    end subroutine calculate_relaxation_mat

    subroutine calculate_polarisation_rotor(time_num, time_num_prop, density_mat, propagator, sizeH, sizeL, &
                pol_i_z, pol_s_z)

        ! Calculate polarisation stroboscopically across multiple rotor periods
        ! Iterative process so number of OMP threads must equal 1

        use omp_lib
        use functions
        implicit none

        integer, intent(in) :: time_num_prop, time_num, sizeH, sizeL
        complex(kind=8), dimension(time_num, sizeL, sizeL), intent(in) :: propagator
        real(kind=8), dimension(sizeH, sizeH), intent(in) :: density_mat
        real(kind=8), dimension(time_num_prop), intent(out) :: pol_i_z, pol_s_z

        integer :: count
        complex(kind=8), dimension(2, 2) :: spin_z, identity_size2

        complex(kind=8), dimension(sizeL, 1) :: density_mat_liouville, temp
        complex(kind=8), dimension(sizeL, sizeL) :: identity_size16, propagator_strobe
        complex(kind=8), dimension(sizeH, sizeH) :: spin2_s_z, spin2_i_z, density_mat_time

        ! Calculate matrices specific to polarisation calculation
        identity_size2 = transpose(reshape((/ 1.D0, 0.D0, 0.D0, 1.D0/), shape(identity_size2)))
        identity_size16 = kron_complex(kron_complex(identity_size2, identity_size2), &
                          kron_complex(identity_size2, identity_size2))

        ! Calculate matrices specific to polarisation calculation
        spin_z = 0.5D0 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))
        spin2_s_z = kron_complex(spin_z, identity_size2)
        spin2_i_z = kron_complex(identity_size2, spin_z)

        density_mat_time = density_mat
        propagator_strobe = identity_size16

        ! Calculate stroboscopic propagator (product of all operators within rotor period)
        do count = 1, time_num
            propagator_strobe = matmul(propagator_strobe, propagator(count, :, :))
        end do

        ! Propogate density matrix using stroboscopic propagator
        do count = 1, time_num_prop

            ! Calculate electronic and nuclear polarisation
            pol_i_z(count) = real(trace_complex(matmul(density_mat_time, spin2_i_z)))
            pol_s_z(count) = real(trace_complex(matmul(density_mat_time, spin2_s_z)))

            ! Transform density matrix (2^N x 2^N to 4^N x 1)
            density_mat_liouville = reshape(density_mat_time, (/sizeL, 1/))

            ! Propagate density matrix
            temp = matmul(propagator_strobe, density_mat_liouville)

            ! Transform density matrix (4^N x 1 to 2^N x 2^N)
            density_mat_time = reshape(temp, (/sizeH, sizeH/))

        end do

    end subroutine calculate_polarisation_rotor

    subroutine calculate_polarisation_sub_rotor(time_num, density_mat, propagator, sizeH, sizeL, &
                pol_i_z_rot, pol_s_z_rot)

        ! Calculate sub rotor polarisation
        ! Iterative process so number of OMP threads must equal 1

        use omp_lib
        use functions
        implicit none

        integer, intent(in) :: time_num, sizeH, sizeL
        complex(kind=8), dimension(time_num, sizeL, sizeL), intent(in) :: propagator
        real(kind=8), dimension(sizeH, sizeH), intent(in) :: density_mat
        real(kind=8), dimension(time_num), intent(out) :: pol_i_z_rot, pol_s_z_rot

        integer :: count
        complex(kind=8), dimension(2, 2) :: spin_z, identity_size2

        complex(kind=8), dimension(sizeL, 1) :: density_mat_liouville, temp
        complex(kind=8), dimension(sizeH, sizeH) :: spin2_s_z, spin2_i_z, density_mat_time

        ! Calculate matrices specific to polarisation calculation
        identity_size2 = transpose(reshape((/ 1.D0, 0.D0, 0.D0, 1.D0/), shape(identity_size2)))
        spin_z = 0.5D0 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))
        spin2_s_z = kron_complex(spin_z, identity_size2)
        spin2_i_z = kron_complex(identity_size2, spin_z)

        density_mat_time = density_mat

        ! Propogate density matrix using stroboscopic propagator
        do count = 1, time_num

            ! Calculate electronic and nuclear polarisation
            pol_i_z_rot(count) = real(trace_complex(matmul(density_mat_time, spin2_i_z)))
            pol_s_z_rot(count) = real(trace_complex(matmul(density_mat_time, spin2_s_z)))

            ! Transform density matrix (2^N x 2^N to 4^N x 1)
            density_mat_liouville = reshape(density_mat_time, (/sizeL, 1/))

            ! Propagate density matrix
            temp = matmul(propagator(count, :, :), density_mat_liouville)

            ! Transform density matrix (4^N x 1 to 2^N x 2^N)
            density_mat_time = reshape(temp, (/sizeH, sizeH/))

        end do

    end subroutine calculate_polarisation_sub_rotor

end module
