module example_subroutine

contains

    subroutine calculate_hamiltonian(hamiltonian)

        use omp_lib
        implicit none

        real(kind = 8), parameter :: PI = 4.D0 * DATAN(1.D0)
        complex(kind = 8), dimension(100, 4, 4), intent(out) :: hamiltonian

        integer :: count
        complex(kind = 8), dimension(4, 4) :: spin2_s_x, spin2_s_z
        complex(kind = 8), dimension(4, 4) :: spin2_i_x, spin2_i_z
        complex(kind = 8), dimension(2, 2) :: spin_x, spin_z, identity_spin1

        ! Identity matrix
        identity_spin1 = transpose(reshape((/ 1.D0, 0.D0, 0.D0, 1.D0/), shape(identity_spin1)))

        ! Pauli matrices
        spin_x = 0.5 * (reshape((/ 0.D0, 1.D0, 1.D0, 0.D0/), shape(spin_x), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))

        ! 4x4 matrices for S operator
        spin2_s_x = kron_complex(spin_x, identity_spin1)
        spin2_s_z = kron_complex(spin_z, identity_spin1)

        ! 4x4 matrices for I operator
        spin2_i_x = kron_complex(identity_spin1, spin_x)
        spin2_i_z = kron_complex(identity_spin1, spin_z)

        do count = 1, 100

            hamiltonian(count, :, :) =  count * 2E6 * spin2_s_z + &
                                        1E6 * spin2_i_z + &
                                        2.D0 * 1E4 * MATMUL(spin2_i_x, spin2_s_z)

        end do

    end subroutine calculate_hamiltonian

    function kron_complex(A, B) result(AB)

        complex(kind = 8), intent(in) :: A(:, :), B(:, :)
        complex(kind = 8) :: AB(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2) )
        INTEGER R, RA, RB, C, CA, CB, I, J

        R = 0
        RA = UBOUND(A, DIM = 1)
        CA = UBOUND(A, DIM = 2)
        RB = UBOUND(B, DIM = 1)
        CB = UBOUND(B, DIM = 2)

        DO I = 1, RA
            C = 0
            DO J = 1, CA
                AB(R + 1 : R + RB, C + 1 : C + CB) = A(I, J) * B
                C = C + CB
            END DO
            R = R + RB
        END DO

    END function kron_complex

end module