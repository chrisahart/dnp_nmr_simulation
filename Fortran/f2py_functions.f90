module f2py_functions

contains

    function trace_real(A) result(C)

        real(kind = 8), dimension (:, :), intent(in) :: A
        real(kind = 8) :: C
        integer :: i
        C = 0

        do i = 1, size(A, 1)
            C = C + A(i, i)
        end do

    end function trace_real

    function trace_complex(A) result(C)

        complex(kind=8), dimension (:, :), intent(in) :: A
        complex(kind=8) :: C
        integer :: i
        C = 0

        do i = 1, size(A, 1)
            C = C + A(i, i)
        end do

    end function trace_complex

    function kron_real(A, B) result(AB)
        ! https://rosettacode.org/wiki/Kronecker_product#Fortran

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

    function expm_complex(t, H) result(expH) ! todo consider reducing precision of calculation

        complex(kind=8), intent(in) :: t
        complex(kind=8), dimension(:, :), intent(in) :: H
        complex(kind=8), dimension(size(H, 1), size(H, 2)) :: expH

        integer, parameter :: ideg = 6
        complex(kind=8), dimension(4 * size(H, 1) * size(H, 2) + ideg + 1) :: wsp
        integer, dimension(size(H, 1)) :: iwsp
        integer :: iexp, ns, iflag, n

        n = size(H, 1)
        call ZGPADM(ideg, n, t, H, n, wsp, size(wsp, 1), iwsp, iexp, ns, iflag)
        expH = reshape(wsp(iexp : iexp + n * n - 1), (/ n, n/))

    end function expm_complex

    function inverse_complex(eig_vector) result(eig_vector_inv)

        integer :: count
        complex(kind=8), dimension(10000, 4, 4), intent(in) :: eig_vector
        complex(kind=8), dimension(10000, 4, 4) :: eig_vector_inv
        complex(kind=8), dimension(4, 4) :: temp

        !integer :: sz1 = 1000, sz2 = 4, sz3 = 8
        integer :: ipiv(4), info_inv
        complex(kind=8) :: work_inv(8)

        eig_vector_inv = 0

        ! Calculate inverse eigenvectors using LAPACK (OpenMP doesn't help here)
        do count = 1, size(eig_vector, 1)
            temp = eig_vector(count, :, :)
            call ZGETRF(4, 4, temp, 4, ipiv, info_inv)
            call ZGETRI(4, temp, 4, ipiv, work_inv, 8, info_inv)
!            call ZGETRF(size(eig_vector, 1), size(eig_vector, 1), temp, size(eig_vector, 1), ipiv, info_inv)
!            call ZGETRI(size(eig_vector, 1), temp, size(eig_vector, 1), ipiv, work_inv, &
!                    size(eig_vector, 1)**2, info_inv)
            eig_vector_inv(count, :, :) = temp
        end do

!        write(6, *) eig_vector_inv(1, 1, 1)

    end function inverse_complex

    subroutine eig_complex(hamiltonian_complex, eigval, eig_vector_complex)

        implicit none

        complex(kind=8), dimension(10000, 4, 4), intent(in) :: hamiltonian_complex
        complex(kind=8), dimension(10000, 4), intent(out) :: eigval
        complex(kind=8), dimension(10000, 4, 4), intent(out) :: eig_vector_complex

        complex(kind=8) :: dummy(4,4), work(8)
        real(kind=8) :: Rwork
        integer :: info, count

        complex(kind=8), dimension(4, 4) :: test1, temp2
        complex(kind=8), dimension(4) :: temp1

        ! Calculate energy, eigenvalues and eigenvectors of intrinsic Hamiltonian
        do count = 1, size(hamiltonian_complex, 1)
            test1 = hamiltonian_complex(count, :, :)
            call ZGEEV('N', 'V', 4, test1, 4, temp1, dummy, 4, temp2, 4, work, 8, Rwork, info)
            eigval(count, :) = temp1
            eig_vector_complex(count, :, :) = temp2
        end do

        write(6,*) eigval(1, :)
        write(6,*) eig_vector_complex(1, :, :)


    end subroutine

end module f2py_functions