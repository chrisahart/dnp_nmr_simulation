module f2py_dynamics

contains

    subroutine calculate_hamiltonian(time_num, time_step, freq_rotor, gtensor, hyperfine_coupling, hyperfine_angles, &
                                    orientation_se, electron_frequency, microwave_frequency, nuclear_frequency, &
                                    hamiltonian)

        use iso_fortran_env
        use omp_lib
        implicit none

        integer, parameter :: WP = REAL64
        complex(kind=8), parameter :: i = (0, 1)
        real(kind=8), parameter :: PI = 4.D0 * DATAN(1.D0)

        integer, intent(in):: time_num
        real(kind=8), dimension(3), intent(in) :: gtensor, hyperfine_angles, orientation_se
        real(kind=8), intent(in) :: hyperfine_coupling, electron_frequency, nuclear_frequency, microwave_frequency
        real(kind=8), intent(in) :: time_step, freq_rotor
        complex(kind=8), intent(out) :: hamiltonian(time_num, 4, 4)

        integer :: count
        complex(kind=8), dimension(4, 4) :: spin2_s_x, spin2_s_y, spin2_s_z, hyperfine_total
        complex(kind=8), dimension(4, 4) :: spin2_i_x, spin2_i_y, spin2_i_z, spin2_identity
        complex(kind=8), dimension(2, 2) :: spin_x, spin_y, spin_z, spin_identity
        complex(kind=8) :: hyperfine_zx, hyperfine_zz, ganisotropy
        real(kind=8) :: gx, gy, gz, ca, cb, cg, sa, sb, sg, r11, r12, r13, r21, r22, r23, r31, r32, r33
        real(kind=8) :: c0, c1, c2, c3, c4

        ! Identity matrix
        spin_identity = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(spin_identity)))

        ! Pauli matrices
        spin_x = 0.5 * (reshape((/ 0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_x), order = (/2, 1/)))
        spin_y = 0.5 * i * (reshape((/ 0.0_WP, -1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_y), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron(spin_x, spin_identity)
        spin2_s_y = kron(spin_y, spin_identity)
        spin2_s_z = kron(spin_z, spin_identity)

        ! 4x4 matrices for I operator
        spin2_i_x = kron(spin_identity, spin_x)
        spin2_i_y = kron(spin_identity, spin_y)
        spin2_i_z = kron(spin_identity, spin_z)

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
