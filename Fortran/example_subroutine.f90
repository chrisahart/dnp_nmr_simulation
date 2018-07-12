module example_subroutine

contains

    subroutine calculate_hamiltonian(hamiltonian)

        use iso_fortran_env
        implicit none

        real(kind = real64), dimension(4, 4), intent(out) :: hamiltonian
        real(kind = real64), dimension(4, 4) :: spin2_s_x, spin2_s_z
        real(kind = real64), dimension(4, 4) :: spin2_i_x, spin2_i_z
        real(kind = real64), dimension(2, 2) :: spin_x, spin_z, identity_spin1

        spin_x = 0.5 * (reshape((/ 0.D0, 1.D0, 1.D0, 0.D0/), shape(spin_x), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))
        identity_spin1 = transpose(reshape((/ 1.D0, 0.D0, 0.D0, 1.D0/), shape(identity_spin1)))

        spin2_s_x = kron_real(spin_x, identity_spin1)
        spin2_s_z = kron_real(spin_z, identity_spin1)

        spin2_i_x = kron_real(identity_spin1, spin_x)
        spin2_i_z = kron_real(identity_spin1, spin_z)

        hamiltonian = 2E6 * spin2_s_z + 1E6 * spin2_i_z + 1E6 * matmul(spin2_i_x, spin2_s_z)

    end subroutine calculate_hamiltonian

    subroutine eig(hamiltonian, eigval, eig_vector)

        use iso_fortran_env
        implicit none

        real(kind = real64), dimension(4, 4), intent(in) :: hamiltonian
        real(kind = real64), intent(out) :: eigval(4), eig_vector(4, 4)

        real(kind = real64) :: dummy1(4), dummy2(4, 4), work(16)
        integer :: info

        call DGEEV('N', 'V', 4, hamiltonian, 4, eigval, dummy1, dummy2, 4, eig_vector, 4, work, 16, info)

    end subroutine eig

    subroutine inv(eig_vector, eig_vector_inv)

        use iso_fortran_env
        implicit none

        real(kind = real64), dimension(4, 4), intent(in) :: eig_vector
        real(kind = real64), dimension(4, 4), intent(out) :: eig_vector_inv

        integer :: info, ipiv(4)
        real(kind=real64) :: work(8)

        eig_vector_inv = eig_vector

        call DGETRF(4, 4, eig_vector_inv, 4, ipiv, info)
        call DGETRI(4, eig_vector_inv, 4, ipiv, work, 4, info)

    end subroutine inv

    function kron_real(A, B) result(AB)

        real(kind = 8), intent(in) :: A(:, :), B(:, :)
        real(kind = 8) :: AB(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2) )
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

    END function kron_real

end module