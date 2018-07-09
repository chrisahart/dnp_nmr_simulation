module f2py_dynamics

contains

    subroutine main(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, &
            orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, microwave_amplitude, &
            t1_nuc,  t1_elec, t2_nuc, t2_elec, temperature, time_num_prop, &
            pol_i_z, pol_s_z, pol_i_z_rot, pol_s_z_rot, energies)

        use iso_fortran_env
        use omp_lib
        use f2py_functions
        implicit none

        integer, intent(in):: time_num, time_num_prop
        real(kind=8), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_se
        real(kind=8), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(kind=8), intent(in) :: time_step, freq_rotor
        real(kind=8), intent(in) :: microwave_amplitude
        real(kind=8), intent(in) :: t1_nuc, t1_elec, t2_nuc, t2_elec, temperature

        complex(kind = 8), dimension(time_num, 16, 16)  :: propagator
        complex(kind = 8), dimension(4, 4) :: density_mat
        complex(kind=8) :: eig_vector(time_num, 4,4), eig_vector_inv(time_num, 4,4)

        real(kind = 8), dimension(time_num_prop), intent(out) :: pol_i_z, pol_s_z
        real(kind = 8), dimension(time_num), intent(out) :: pol_i_z_rot, pol_s_z_rot
        real(kind = 8), dimension(time_num, 4), intent(out) :: energies

        call calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, &
            hyperfine_coupling, hyperfine_angles, &
                                    orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, &
            energies, eig_vector, eig_vector_inv, density_mat)

        call liouville_propagator(time_num, time_step, electron_frequency, nuclear_frequency, microwave_amplitude, &
            t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, eig_vector, eig_vector_inv, energies, propagator)

        call calculate_polarisation_rotor(time_num, time_num_prop, density_mat, propagator, pol_i_z, pol_s_z)

        call calculate_polarisation_sub_rotor(time_num, density_mat, propagator, pol_i_z_rot, pol_s_z_rot)

    end subroutine main

    subroutine calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, &
            hyperfine_coupling, hyperfine_angles, &
                                    orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, &
            energies, eig_vector, eig_vector_inv, density_mat)

        use iso_fortran_env
        use omp_lib
        use f2py_functions
        implicit none

        integer, parameter :: WP = REAL64
        complex(kind=8), parameter :: i = (0, 1)
        real(kind=8), parameter :: PI = 4.D0 * DATAN(1.D0), Planck = 6.62607004E-34, Boltzmann = 1.38064852E-23

        integer, intent(in):: time_num
        real(kind=8), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_se
        real(kind=8), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(kind=8), intent(in) :: time_step, freq_rotor
        real(kind=8), intent(out) :: energies(time_num, 4)
        complex(kind=8), intent(out) :: eig_vector(time_num, 4,4), eig_vector_inv(time_num, 4,4), density_mat(4, 4)
        complex(kind=8) :: temp(time_num, 4,4)

        integer :: count, count2
        complex(kind=8), dimension(4) :: test3
        complex(kind=8), dimension(4, 4) :: spin2_s_x, spin2_s_y, spin2_s_z, hyperfine_total, test1, test2
        complex(kind=8), dimension(4, 4) :: spin2_i_x, spin2_i_y, spin2_i_z, spin2_identity
        complex(kind=8), dimension(2, 2) :: spin_x, spin_y, spin_z, identity_spin1
        complex(kind=8) :: hyperfine_zx, hyperfine_zz, ganisotropy
        real(kind=8) :: gx, gy, gz, ca, cb, cg, sa, sb, sg, r11, r12, r13, r21, r22, r23, r31, r32, r33
        real(kind=8) :: c0, c1, c2, c3, c4

        real(kind=8) :: Rwork
        complex(kind=8) :: eigval(time_num, 4), dummy(1,1), work(8), work_inv(8)
        integer :: info, info_inv
        complex(kind=8), dimension(time_num, 4,4) :: hamiltonian
        double precision wtime

        complex(kind=8) :: N(6,6),VL(1,1),W(6),WORK2(12),VR(6,6),NN(6,6)
        real(kind=8) :: RWORK2(12)
        integer :: INFO2
        integer :: ipiv(4)

        complex(kind=8), dimension(4, 4) :: hamiltonian_ideal, boltzmann_factors_mat
        complex(kind=8), dimension(4) :: boltzmann_factors, energies_ideal

        ! Identity matrix
        identity_spin1 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_spin1)))

        ! Pauli matrices
        spin_x = 0.5 * (reshape((/ 0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_x), order = (/2, 1/)))
        spin_y = 0.5 * i * (reshape((/ 0.0_WP, -1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_y), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron_complex(spin_x, identity_spin1)
        spin2_s_y = kron_complex(spin_y, identity_spin1)
        spin2_s_z = kron_complex(spin_z, identity_spin1)

        ! 4x4 matrices for I operator
        spin2_i_x = kron_complex(identity_spin1, spin_x)
        spin2_i_y = kron_complex(identity_spin1, spin_y)
        spin2_i_z = kron_complex(identity_spin1, spin_z)

        ! Calculate g-anisotropy
        gx = electron_frequency * gtensor(1)
        gy = electron_frequency * gtensor(2)
        gz = electron_frequency * gtensor(3)

        ca = cos(orientation_se(1))
        cb = cos(orientation_se(2))
        cg = cos(orientation_se(3))
        sa = sin(orientation_se(1))
        sb = sin(orientation_se(2))
        sg = sin(orientation_se(3))

        r11 = ca * cb * cg - sa * sg
        r12 = sa * cb * cg + ca * sg
        r13 = -sb * cg
        r21 = -ca * cb * sg - sa * cg
        r22 = -sa * cb * sg + ca * cg
        r23 = sb * sg
        r31 = ca * sb
        r32 = sa * sb
        r33 = cb

        c0 = 1.0_WP / 3.0_WP * (gx + gy + gz)
        c1 = 2.0_WP * 2.0_WP ** (1.0_WP/2.0_WP) / 3.0_WP * (gx * r11 * r31 + gy * r12 * r32 + gz * r13 * r33)
        c2 = 2.0_WP ** (1.0_WP/2.0_WP) / 3.0_WP * (gx * r21 * r31 + gy * r22 * r32 + gz * r23 * r33)
        c3 = 1.0_WP / 3.0_WP * (gx * (r11 ** 2.0_WP - r21 ** 2.0_WP) + gy * (r12 ** 2.0_WP - r22 ** 2.0_WP) + &
                gz * (r13 ** 2 - r23 ** 2))
        c4 = 2.0_WP / 3.0_WP * (gx * r11 * r21 + gy * r22 * r12 + gz * r13 * r23)

        !!$omp parallel do default(private) &
        !!$omp& shared(c0, c1, c2, c3, c4, freq_rotor, time_step, hyperfine_angles, hyperfine_coupling) &
        !!$omp& shared(hamiltonian, microwave_frequency, nuclear_frequency) &
        !!$omp& shared(spin2_i_x, spin2_i_y, spin2_i_z, spin2_s_z)
        do count = 1, time_num

            ! Calculate time dependent hyperfine
            hyperfine_zz = hyperfine_coupling * (-0.5_WP * (sin(hyperfine_angles(2)) ** 2.0_WP) * &
                    cos(2.0_WP * (2.0_WP * PI * freq_rotor * (count - 1) * time_step + &
                            hyperfine_angles(3))) + 2.0_WP ** 2.0_WP * sin(hyperfine_angles(2)) * &
                    cos(hyperfine_angles(2)) * &
                    cos(2 * PI * freq_rotor * (count - 1) * time_step + hyperfine_angles(3)))

            hyperfine_zx = hyperfine_coupling * (-0.5_WP * sin(hyperfine_angles(2)) * cos(hyperfine_angles(2)) * &
                    cos(2_WP * PI * freq_rotor * (count - 1) * time_step + &
                            hyperfine_angles(3)) - (2.0_WP ** 0.5_WP / 4.0_WP) * (sin(hyperfine_angles(2)) ** 2.0_WP)*&
                    cos(2.0_WP * (2.0_WP * PI * freq_rotor * (count - 1) * time_step + hyperfine_angles(3))) + &
                    (2.0_WP ** 0.5_WP / 4.0_WP) * (3.0_WP * (cos(hyperfine_angles(2)) ** 2.0_WP) - 1.0_WP))

            hyperfine_total = 2.0_WP * hyperfine_zx * MATMUL(spin2_i_x, spin2_s_z) + &
                    2.0_WP * hyperfine_zz * MATMUL(spin2_i_z, spin2_s_z)

            ! Calculate time dependent g-anisotropy frequency
            ganisotropy = c0 + c1 * cos(2.0_WP * PI * freq_rotor * (count - 1) * time_step) + &
                    c2 * sin(2.0_WP * PI * freq_rotor * (count - 1) * time_step) + &
                    c3 * cos(4.0_WP * PI * freq_rotor * (count - 1) * time_step) + &
                    c4 * sin(4.0_WP * PI * freq_rotor * (count - 1) * time_step)


            ! Calculate time dependent hamiltonian
            hamiltonian(count, :, :) = (ganisotropy - microwave_frequency) * spin2_s_z + &
                    nuclear_frequency * spin2_i_z + &
                    hyperfine_total

        end do
        !!$omp end parallel do

        ! Calculate eigenvalues and eigenvectors using LAPACK (OpenMP doesn't help here)
        do count2 = 1, time_num
            test1 = hamiltonian(count2, :, :)
            call ZGEEV('N', 'V', 4, test1, 4, eigval(count2, :), dummy, 1,  &
                    eig_vector(count2, :, :), 4, work, 8, Rwork, info)
            write(6, *) eig_vector(1, 1, 1)
        end do
        energies = real(eigval)
        write(6, *) eig_vector(1, 1, 1)

        ! Getting some very odd precision issues, eigvector=1.0000000000000000 if inverse eigenvector loop not present
        ! If present then eigvector=0.99999650020944664
        ! Convert complex to real, then look into precision

         ! Calculate inverse eigenvectors using LAPACK (OpenMP doesn't help here)
        do count2 = 1, size(eig_vector, 1)
            test2 = eig_vector(count2, :, :)
            call ZGETRF(4, 4, test2, 4, ipiv, info_inv)
            call ZGETRI(4, test2, 4, ipiv, work_inv, 8, info_inv)
            eig_vector_inv(count2, :, :) = test2
        end do
!        write(6, *) eig_vector(1, 1, 1)
!        write(6, *) eig_vector_inv(1, 1, 1)
!
!!        temp = inverse_complex(eig_vector)
!        write(6, *) temp(1, 1, 1)
!        eig_vector_inv = temp
!        write(6, *) eig_vector_inv(1, 1, 1)

        ! Calculate thermal density matrix from idealised Hamiltonian
        hamiltonian_ideal = electron_frequency * spin2_s_z + nuclear_frequency * spin2_i_z
        energies_ideal(1) = hamiltonian_ideal(1, 1)
        energies_ideal(2) = hamiltonian_ideal(2, 2)
        energies_ideal(3) = hamiltonian_ideal(3, 3)
        energies_ideal(4) = hamiltonian_ideal(4, 4)

        ! Calculate initial Zeeman basis density matrix from Boltzmann factors
        do count = 1, 4
            boltzmann_factors(count) = exp(-(Planck * energies_ideal(count)) / (Boltzmann * 100))
        end do
        boltzmann_factors_mat(1, 1) = boltzmann_factors(1)
        boltzmann_factors_mat(2, 2) = boltzmann_factors(2)
        boltzmann_factors_mat(3, 3) = boltzmann_factors(3)
        boltzmann_factors_mat(4, 4) = boltzmann_factors(4)
        density_mat = 0
        density_mat = (1/sum(boltzmann_factors)) * boltzmann_factors_mat

    end subroutine calculate_hamiltonian

    subroutine liouville_propagator(time_num, time_step, electron_frequency, freq_nuclear_1, microwave_amplitude, &
            t1_nuc, t1_elec, t2_nuc, t2_elec, temperature, eigvectors, eigvectors_inv, energies, propagator)

        use iso_fortran_env
        use omp_lib
        use f2py_functions
        implicit none

        integer, parameter :: WP = REAL64
        complex(kind=8), parameter :: i = (0, 1)
        real(kind=8), parameter :: PI = 4.D0 * DATAN(1.D0), Planck = 6.62607004E-34, Boltzmann = 1.38064852E-23

        integer, intent(in):: time_num
        complex(kind=8), dimension(:, :, :), intent(in) :: eigvectors, eigvectors_inv
        real(kind=8), dimension(:, :), intent(in) :: energies
        real(kind=8), intent(in) :: electron_frequency, freq_nuclear_1, microwave_amplitude, time_step
        real(kind=8), intent(in) :: t1_nuc, t1_elec, t2_nuc, t2_elec, temperature
        complex(kind=8), intent(out) :: propagator(time_num, 16, 16)

        integer :: count
        complex(kind=8), dimension(16, 16) :: identity_size16, hamiltonian_liouville, relax_mat, mat_exp
        complex(kind=8), dimension(16, 16) :: eigvectors_liouville, eigvectors_inv_liouville, liouvillian, exponent
        complex(kind=8), dimension(16, 16) :: spin2_i_p_tl, spin2_i_m_tl, spin2_s_p_tl, spin2_s_m_tl
        complex(kind=8), dimension(16, 16) :: relax_t2_elec, relax_t2_nuc, relax_t1
        complex(kind=8), dimension(4, 4) :: spin2_s_z_t, spin2_s_p_t, spin2_s_m_t, spin2_i_z_t, spin2_i_p_t, spin2_i_m_t
        complex(kind=8), dimension(4, 4) :: spin2_s_x, spin2_s_y, spin2_s_z, total_hamiltonian
        complex(kind=8), dimension(4, 4) :: spin2_s_p, spin2_s_m, spin2_i_p, spin2_i_m, test1
        complex(kind=8), dimension(4, 4) :: microwave_hamiltonian_init, microwave_hamiltonian, energy_mat
        complex(kind=8), dimension(4, 4) :: spin2_i_x, spin2_i_y, spin2_i_z, identity_size4
        complex(kind=8), dimension(2, 2) :: spin_x, spin_y, spin_z, identity_size2
        real(kind=8) :: p_e, p_n, gnp, gnm, gep, gem
        complex(kind=8) :: timescale

        ! Identity matrix
        identity_size2 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_size2)))
        identity_size4 = kron_complex(identity_size2, identity_size2)
        identity_size16 = kron_complex(kron_complex(identity_size2, identity_size2), &
                kron_complex(identity_size2, identity_size2))

        ! Pauli matrices
        spin_x = 0.5 * (reshape((/ 0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_x), order = (/2, 1/)))
        spin_y = 0.5 * i * (reshape((/ 0.0_WP, -1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_y), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron_complex(spin_x, identity_size2)
        spin2_s_y = kron_complex(spin_y, identity_size2)
        spin2_s_z = kron_complex(spin_z, identity_size2)
        spin2_s_p = spin2_s_x + i * spin2_s_y
        spin2_s_m = spin2_s_x - i * spin2_s_y

        ! 4x4 matrices for I operator
        spin2_i_x = kron_complex(identity_size2, spin_x)
        spin2_i_y = kron_complex(identity_size2, spin_y)
        spin2_i_z = kron_complex(identity_size2, spin_z)
        spin2_i_p = spin2_i_x + i * spin2_i_y
        spin2_i_m = spin2_i_x - i * spin2_i_y

        ! Calculate variables for Liouville space relaxation
        p_e = tanh(0.5 * electron_frequency * (Planck / (Boltzmann * temperature)))
        p_n = tanh(0.5 * freq_nuclear_1 * (Planck / (Boltzmann * temperature)))
        gnp = 0.5 * (1 - p_n) * (1 / (1 * t1_nuc))
        gnm = 0.5 * (1 + p_n) * (1 / (1 * t1_nuc))
        gep = 0.5 * (1 - p_e) * (1 / (1 * t1_elec))
        gem = 0.5 * (1 + p_e) * (1 / (1 * t1_elec))

        ! Calculate initial microwave Hamiltonian
        microwave_hamiltonian_init = microwave_amplitude * spin2_s_x

        !!$omp parallel do default(private) &
        !!$omp& shared(microwave_hamiltonian_init, eigvectors, eigvectors_inv, energies, time_step)&
        !!$omp& shared(t2_elec, t2_nuc, gnp, gnm, gep, gem, propagator) &
        !!$omp& shared(spin2_s_z, spin2_s_p, spin2_s_m, spin2_i_z, spin2_i_p, spin2_i_m, identity_size4, identity_size16)
        do count = 1, time_num

            ! Transform microwave Hamiltonian into time dependent basis
            microwave_hamiltonian = matmul(eigvectors_inv(count, :, :), matmul(microwave_hamiltonian_init, &
                    eigvectors(count, :, :)))

            ! Calculate total Hamiltonian
            energy_mat = identity_size4
            energy_mat(1, 1) = energies(count, 1)
            energy_mat(2, 2) = energies(count, 2)
            energy_mat(3, 3) = energies(count, 3)
            energy_mat(4, 4) = energies(count, 4)

            ! Calculate total Hamiltonian
            total_hamiltonian = energy_mat + microwave_hamiltonian

            ! Transform Hilbert space Hamiltonian into Liouville space
            hamiltonian_liouville = kron_complex(total_hamiltonian, identity_size4) - &
                    kron_complex(identity_size4, transpose(total_hamiltonian))

            ! Transform spin matrices into time dependent Hilbert space basis
            spin2_s_z_t = matmul(eigvectors_inv(count, :, :), matmul(spin2_s_z, eigvectors(count, :, :)))
            spin2_s_p_t = matmul(eigvectors_inv(count, :, :), matmul(spin2_s_p, eigvectors(count, :, :)))
            spin2_s_m_t = matmul(eigvectors_inv(count, :, :), matmul(spin2_s_m, eigvectors(count, :, :)))
            spin2_i_z_t = matmul(eigvectors_inv(count, :, :), matmul(spin2_i_z, eigvectors(count, :, :)))
            spin2_i_p_t = matmul(eigvectors_inv(count, :, :), matmul(spin2_i_p, eigvectors(count, :, :)))
            spin2_i_m_t = matmul(eigvectors_inv(count, :, :), matmul(spin2_i_m, eigvectors(count, :, :)))

            ! Transform spin matrices into time dependent Liouville space basis
            spin2_i_p_tl = kron_complex(spin2_i_p_t, transpose(spin2_i_m_t)) - 0.5_WP * identity_size16 + 0.5_WP * ( &
                           kron_complex(spin2_i_z_t, identity_size4) + &
                           kron_complex(identity_size4, transpose(spin2_i_z_t)))

            spin2_i_m_tl = kron_complex(spin2_i_m_t, transpose(spin2_i_p_t)) - 0.5_WP * identity_size16 - 0.5_WP * ( &
                            kron_complex(spin2_i_z_t, identity_size4) + &
                            kron_complex(identity_size4, transpose(spin2_i_z_t)))

            spin2_s_p_tl = kron_complex(spin2_s_p_t, transpose(spin2_s_m_t)) - 0.5_WP * identity_size16 + 0.5_WP * ( &
                           kron_complex(spin2_s_z_t, identity_size4) + &
                           kron_complex(identity_size4, transpose(spin2_s_z_t)))

            spin2_s_m_tl = kron_complex(spin2_s_m_t, transpose(spin2_s_p_t)) - 0.5_WP * identity_size16 - 0.5_WP * ( &
                           kron_complex(spin2_s_z_t, identity_size4) + &
                           kron_complex(identity_size4, transpose(spin2_s_z_t)))

            ! Calculate time dependent Liouville space relaxation matrix
            relax_t2_elec = (1.0_WP / t2_elec) * (kron_complex(spin2_s_z_t, transpose(spin2_s_z_t)) - &
                    0.5_WP * 0.5_WP * identity_size16)
            relax_t2_nuc = (1.0_WP / t2_nuc) * (kron_complex(spin2_i_z_t, transpose(spin2_i_z_t)) - &
                    0.5_WP * 0.5_WP * identity_size16)
            relax_t1 = gep * spin2_s_p_tl + gem * spin2_s_m_tl + gnp * spin2_i_p_tl + gnm * spin2_i_m_tl
            relax_mat = relax_t2_elec + relax_t2_nuc + relax_t1

            ! Calculate Liouville space eigenvectors
            eigvectors_liouville = kron_complex(eigvectors(count, :, :), eigvectors(count, :, :))
            eigvectors_inv_liouville = kron_complex(eigvectors_inv(count, :, :), eigvectors_inv(count, :, :))

            ! Calculate Liouville space propagator
            liouvillian = hamiltonian_liouville + i * relax_mat
            exponent = -i * liouvillian * time_step
            timescale = 1.0_WP
            mat_exp = expm_complex(timescale, exponent)
            propagator(count, :, :) = matmul(eigvectors_inv_liouville, matmul(mat_exp, eigvectors_liouville))

        end do
        !!$omp end parallel do

        !propagator = real(propagator)

    end subroutine liouville_propagator


    subroutine calculate_polarisation_rotor(time_num, time_num_prop, density_mat, propagator, pol_i_z, pol_s_z)

        use iso_fortran_env
        use omp_lib
        use f2py_functions
        implicit none

        integer, parameter :: WP = REAL64

        integer, intent(in) :: time_num_prop, time_num
        complex(kind = 8), dimension(:, :, :), intent(in) :: propagator
        complex(kind = 8), dimension(:, :), intent(in) :: density_mat
        real(kind = 8), dimension(time_num_prop), intent(out) :: pol_i_z, pol_s_z

        integer :: count
        complex(kind = 8), dimension(16, 1) :: density_mat_liouville, temp
        complex(kind = 8), dimension(16, 16) :: identity_size16, propagator_strobe
        complex(kind = 8), dimension(4, 4) :: spin2_s_z, spin2_i_z, density_mat_time
        complex(kind = 8), dimension(2, 2) :: spin_z, identity_size2

        ! Calculate matrices specific to polarisation calculation
        identity_size2 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_size2)))
        identity_size16 = kron_complex(kron_complex(identity_size2, identity_size2), &
                kron_complex(identity_size2, identity_size2))
        spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))
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
            pol_i_z(count) = trace_complex(matmul(density_mat_time, spin2_i_z))
            pol_s_z(count) = trace_complex(matmul(density_mat_time, spin2_s_z))

            ! Transform density matrix (2^N x 2^N to 4^N x 1)
            density_mat_liouville = reshape(density_mat_time, (/4 ** 2, 1/))

            ! Propagate density matrix
            temp = matmul(propagator_strobe, density_mat_liouville)

            ! Transform density matrix (4^N x 1 to 2^N x 2^N)
            density_mat_time = reshape(temp, (/2 ** 2, 2 ** 2/))

        end do

    end subroutine calculate_polarisation_rotor

    subroutine calculate_polarisation_sub_rotor(time_num, density_mat, propagator, pol_i_z_rot, pol_s_z_rot)

        use iso_fortran_env
        use omp_lib
        use f2py_functions
        implicit none

        integer, parameter :: WP = REAL64

        integer, intent(in) :: time_num
        complex(kind = 8), dimension(:, :, :), intent(in) :: propagator
        complex(kind = 8), dimension(:, :), intent(in) :: density_mat
        real(kind = 8), dimension(time_num), intent(out) :: pol_i_z_rot, pol_s_z_rot

        integer :: count
        complex(kind = 8), dimension(16, 1) :: density_mat_liouville, temp
        complex(kind = 8), dimension(4, 4) :: spin2_s_z, spin2_i_z, density_mat_time
        complex(kind = 8), dimension(2, 2) :: spin_z, identity_size2

        ! Calculate matrices specific to polarisation calculation
        identity_size2 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_size2)))
        spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))
        spin2_s_z = kron_complex(spin_z, identity_size2)
        spin2_i_z = kron_complex(identity_size2, spin_z)

        density_mat_time = density_mat

        do count = 1, time_num

            ! Calculate electronic and nuclear polarisation
            pol_i_z_rot(count) = trace_complex(matmul(density_mat_time, spin2_i_z))
            pol_s_z_rot(count) = trace_complex(matmul(density_mat_time, spin2_s_z))

            ! Transform density matrix (2^N x 2^N to 4^N x 1)
            density_mat_liouville = reshape(density_mat_time, (/4 ** 2, 1/))

            ! Propagate density matrix
            temp = matmul(propagator(count, :, :), density_mat_liouville)

            ! Transform density matrix (4^N x 1 to 2^N x 2^N)
            density_mat_time = reshape(temp, (/2 ** 2, 2 ** 2/))

        end do

    end subroutine calculate_polarisation_sub_rotor

end module
