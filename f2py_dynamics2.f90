module f2py_dynamics

contains

    subroutine calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, hamiltonian_ideal, &
            hyperfine_coupling, hyperfine_angles, &
                                    orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, &
                                    hamiltonian)

        use iso_fortran_env
        use omp_lib
        implicit none

        integer, parameter :: WP = REAL64
        complex(kind=8), parameter :: i = (0, 1)
        real(kind=8), parameter :: PI = 4.D0 * DATAN(1.D0)

        integer, intent(in):: time_num
        real(kind=8), dimension(:, :, :), intent(in) :: hamiltonian_ideal
        real(kind=8), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_se
        real(kind=8), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(kind=8), intent(in) :: time_step, freq_rotor
        complex(kind=8), intent(out) :: hamiltonian(time_num, 4, 4)

        integer :: count
        complex(kind=8), dimension(4, 4) :: spin2_s_x, spin2_s_y, spin2_s_z, hyperfine_total
        complex(kind=8), dimension(4, 4) :: spin2_i_x, spin2_i_y, spin2_i_z, spin2_identity
        complex(kind=8), dimension(2, 2) :: spin_x, spin_y, spin_z, identity_spin1
        complex(kind=8) :: hyperfine_zx, hyperfine_zz, ganisotropy
        real(kind=8) :: gx, gy, gz, ca, cb, cg, sa, sb, sg, r11, r12, r13, r21, r22, r23, r31, r32, r33
        real(kind=8) :: c0, c1, c2, c3, c4

        write(6, *), shape(hamiltonian_ideal)

        ! Identity matrix
        identity_spin1 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_spin1)))

        ! Pauli matrices
        spin_x = 0.5 * (reshape((/ 0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_x), order = (/2, 1/)))
        spin_y = 0.5 * i * (reshape((/ 0.0_WP, -1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_y), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron(spin_x, identity_spin1)
        spin2_s_y = kron(spin_y, identity_spin1)
        spin2_s_z = kron(spin_z, identity_spin1)

        ! 4x4 matrices for I operator
        spin2_i_x = kron(identity_spin1, spin_x)
        spin2_i_y = kron(identity_spin1, spin_y)
        spin2_i_z = kron(identity_spin1, spin_z)

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

        !$omp parallel do default(private) &
        !$omp& shared(c0, c1, c2, c3, c4, freq_rotor, time_step, hyperfine_angles, hyperfine_coupling) &
        !$omp& shared(hamiltonian, microwave_frequency, nuclear_frequency) &
        !$omp& shared(spin2_i_x, spin2_i_y, spin2_i_z, spin2_s_z)
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
        !$omp end parallel do

    end subroutine calculate_hamiltonian

    subroutine liouville_propagator(time_num, time_step, electron_frequency, freq_nuclear_1, microwave_amplitude, &
            t1_nuc, t1_elec, temperature, eigvectors, eigvectors_inv, energies, propagator)

        use iso_fortran_env
        use omp_lib
        implicit none

        integer, parameter :: WP = REAL64
        complex(kind=8), parameter :: i = (0, 1)
        real(kind=8), parameter :: PI = 4.D0 * DATAN(1.D0), Planck = 6.62607004E-34, Boltzmann = 1.38064852E-23

        integer, intent(in):: time_num
        complex(kind=8), dimension(time_num, 4, 4), intent(in) :: eigvectors, eigvectors_inv
        real(kind=8), dimension(time_num, 4), intent(in) :: energies
        real(kind=8), intent(in) :: electron_frequency, freq_nuclear_1, microwave_amplitude, time_step
        real(kind=8), intent(in) :: t1_nuc, t1_elec, temperature
        complex(kind=8), intent(out) :: propagator(time_num, 16, 16)

        integer :: count
        complex(kind=8), dimension(16, 16) :: identity_size16, hamiltonian_liouville, relax_mat
        complex(kind=8), dimension(16, 16) :: eigvectors_liouville, eigvectors_inv_liouville, liouvillian
        complex(kind=8), dimension(4, 4) :: spin2_s_x, spin2_s_y, spin2_s_z, total_hamiltonian
        complex(kind=8), dimension(4, 4) :: microwave_hamiltonian_init, microwave_hamiltonian, energy_mat
        complex(kind=8), dimension(4, 4) :: spin2_i_x, spin2_i_y, spin2_i_z, identity_size4
        complex(kind=8), dimension(2, 2) :: spin_x, spin_y, spin_z, identity_size2
        real(kind=8) :: p_e, p_n, gnp, gnm, gep, gem

        ! Identity matrix
        identity_size2 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_size2)))
        identity_size4 = kron(identity_size2, identity_size2)
        identity_size16 = kron(kron(identity_size2, identity_size2), kron(identity_size2, identity_size2))

        ! Pauli matrices
        spin_x = 0.5 * (reshape((/ 0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_x), order = (/2, 1/)))
        spin_y = 0.5 * i * (reshape((/ 0.0_WP, -1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_y), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron(spin_x, identity_size2)
        spin2_s_y = kron(spin_y, identity_size2)
        spin2_s_z = kron(spin_z, identity_size2)

        ! 4x4 matrices for I operator
        spin2_i_x = kron(identity_size2, spin_x)
        spin2_i_y = kron(identity_size2, spin_y)
        spin2_i_z = kron(identity_size2, spin_z)

        ! Calculate variables for Liouville space relaxation
        p_e = tanh(0.5 * electron_frequency * (Planck / (Boltzmann * temperature)))
        p_n = tanh(0.5 * freq_nuclear_1 * (Planck / (Boltzmann * temperature)))
        gnp = 0.5 * (1 - p_n) * (1 / (1 * t1_nuc))
        gnm = 0.5 * (1 + p_n) * (1 / (1 * t1_nuc))
        gep = 0.5 * (1 - p_e) * (1 / (1 * t1_elec))
        gem = 0.5 * (1 + p_e) * (1 / (1 * t1_elec))

        ! Calculate initial microwave Hamiltonian
        microwave_hamiltonian_init = microwave_amplitude * spin2_s_x

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
            hamiltonian_liouville = kron(total_hamiltonian, identity_size4) - &
                    kron(identity_size4, transpose(total_hamiltonian))

            ! Calculate time dependent Liouville space relaxation matrix
            relax_mat = 0.5 * identity_size16

            ! Calculate Liouville space eigenvectors
            eigvectors_liouville = kron(eigvectors(count, :, :), eigvectors(count, :, :))
            eigvectors_inv_liouville = kron(eigvectors_inv(count, :, :), eigvectors_inv(count, :, :))

            ! Calculate Liouville space propagator
            liouvillian = hamiltonian_liouville + i * relax_mat
            propagator(count, :, :) = expm(eigvectors_inv_liouville, eigvectors_liouville)
!            propagator(count, :, :) = matmul(eigvectors_inv_liouville, &
!                                         matmul(expm(-i * liouvillian * time_step), eigvectors_liouville))
        end do

    end subroutine liouville_propagator

    function expm(A, B) result(C)

        complex(kind=8), dimension (:, :), intent(in) :: A, B
        complex(kind=8), dimension (:, :), allocatable :: C

        allocate(C(size(A, 1), size(A, 2)))
        C = matmul(A, B)

    end function expm

    function kron(A, B) result(C)

        complex(kind=8), dimension (:, :), intent(in) :: A, B
        complex(kind=8), dimension (:, :), allocatable :: C
        integer :: i = 0, j = 0, k = 0, l = 0
        integer :: m = 0, n = 0, p = 0, q = 0

        allocate(C(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2)))
        C = 0

        do i = 1, size(A, 1)
            do j = 1, size(A, 2)
                n = (i - 1) * size(B, 1) + 1
                m = n + size(B, 1)
                p = (j - 1) * size(B, 2) + 1
                q = p + size(B, 2)
                C(n : m, p : q) = A(i, j) * B
            enddo
        enddo

    end function kron

end module
