module functions

    ! This module contains an assortment of common Fortran functions and subroutines for 2D and 3D arrays.

contains

    function trace_real(A) result(C)
        ! Calculate trace of real matrix A

        implicit none

        real(kind = 8), dimension (:, :), intent(in) :: A
        real(kind = 8) :: C
        integer :: i

        C = 0
        do i = 1, size(A, 1)
            C = C + A(i, i)
        end do

    end function trace_real

    function trace_complex(A) result(C)
        ! Calculate trace of complex matrix A

        implicit none

        complex(kind=8), dimension (:, :), intent(in) :: A
        complex(kind=8) :: C
        integer :: i

        C = 0
        do i = 1, size(A, 1)
            C = C + A(i, i)
        end do

    end function trace_complex

    function kron_real(A, B) result(AB)
        ! Calculate Kronecker product of real matrices A and B, adapted from Rosetta Code

        use omp_lib
        implicit none

        real(kind = 8), intent(in) :: A(:, :), B(:, :)
        real(kind = 8) :: AB(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2) )
        integer :: R, RA, RB, C, CA, CB, I, J

        R = 0
        RA = size(A, DIM = 1)
        CA = size(A, DIM = 2)
        RB = size(B, DIM = 1)
        CB = size(B, DIM = 2)

        !!$omp parallel do default(private) &
        !!$omp& shared(A, B, AB, R, CA, RB, CB)
        do I = 1, RA
            C = 0
            do J = 1, CA
                AB(R + 1 : R + RB, C + 1 : C + CB) = A(I, J) * B
                C = C + CB
            end do
            R = R + RB
        end do
        !!$omp end parallel do

    end function kron_real

    function kron_complex(A, B) result(AB)
        ! Calculate Kronecker product of complex matrices A and B, adapted from Rosetta Code

        use omp_lib
        implicit none

        complex(kind = 8), intent(in) :: A(:, :), B(:, :)
        complex(kind = 8) :: AB(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2) )
        integer :: R, RA, RB, C, CA, CB, I, J

        integer threads
        RA = size(A, DIM = 1)

        !!call omp_set_num_threads(1)
        !!$omp parallel do default(private) &
        !!$omp& shared(A, B, AB)
        do I = 1, RA
             R = 0
            RA = size(A, DIM = 1)
            CA = size(A, DIM = 2)
            RB = size(B, DIM = 1)
            CB = size(B, DIM = 2)
            !threads = omp_get_num_threads()
            !write(6,*) 'threads',  threads

            C = 0
            do J = 1, CA
                AB(R + 1 : R + RB, C + 1 : C + CB) = A(I, J) * B
                C = C + CB
            end do
            R = R + RB
        end do
        !!$omp end parallel do

    end function kron_complex

    function inverse_real(A) result(B)
        ! Calculate inverse of 3D real matrix A using LAPACK

        implicit none

        real(kind=8), dimension(:, :, :), intent(in) :: A
        real(kind=8), dimension(size(A, 1), size(A, 2), size(A, 3)) :: B
        real(kind=8), dimension(size(A, 2), size(A, 3)) :: temp
        real(kind=8) :: work(size(A, 2)*2)

        integer :: count
        integer :: ipiv(4), info

        do count = 1, size(A, 1)
            temp = A(count, :, :)
            call DGETRF(size(A, 2), size(A, 2), temp, size(A, 2), ipiv, info)
            call DGETRI(size(A, 2), temp, size(A, 2), ipiv, work, size(A, 1)*2, info)
            B(count, :, :) = temp
        end do

    end function inverse_real

    function inverse_complex(A) result(B)
        ! Calculate inverse of 3D complex matrix A using LAPACK

        implicit none

        complex(kind=8), dimension(:, :, :), intent(in) :: A
        complex(kind=8), dimension(size(A, 1), size(A, 2), size(A, 3)) :: B
        complex(kind=8), dimension(size(A, 2), size(A, 3)) :: temp
        complex(kind=8) :: work(size(A, 2)*2)

        integer :: count
        integer :: ipiv(4), info

        do count = 1, size(A, 1)
            temp = A(count, :, :)
            call ZGETRF(size(A, 2), size(A, 2), temp, size(A, 2), ipiv, info)
            call ZGETRI(size(A, 2), temp, size(A, 2), ipiv, work, size(A, 1)*2, info)
            B(count, :, :) = temp
        end do

    end function inverse_complex

    function expm_complex(A) result(B)
        ! Calculate matrix exponential of complex matrix A using Expokit

        implicit none

        complex(kind=8), dimension(:, :), intent(in) :: A
        complex(kind=8), dimension(size(A, 1), size(A, 2)) :: B

        integer, parameter :: ideg = 2 ! Pade approximation, 6 is reccomended but 2 appears to be stable
        complex(kind=8) :: t = 1.D0
        complex(kind=8), dimension(4 * size(A, 1) * size(A, 2) + ideg + 1) :: wsp
        integer, dimension(size(A, 1)) :: iwsp
        integer :: iexp, ns, iflag

        call ZGPADM(ideg, size(A, 1), t, A, size(A, 1), wsp, size(wsp, 1), iwsp, iexp, ns, iflag)
        B = reshape(wsp(iexp : iexp + size(A, 1) * size(A, 1) - 1), (/ size(A, 1), size(A, 2)/))

    end function expm_complex

    subroutine eig_real(A, eig_val, eig_vect)
        ! Calculate eigenvalues and eigenvectors of 3D real matrix A using LAPACK

        implicit none

        real(kind=8), dimension(:, :, :), intent(in) :: A
        real(kind=8), dimension(size(A, 1), size(A, 2)), intent(out) :: eig_val
        real(kind=8), dimension(size(A, 1), size(A, 2), size(A, 3)), intent(out) :: eig_vect

        real(kind=8) :: dummy_vect(size(A, 2), size(A, 2)), dummy_val(size(A, 2)), work(size(A, 2)*4)
        integer :: info, count

        real(kind=8), dimension(size(A, 2), size(A, 2)) :: test1, temp2
        real(kind=8), dimension(size(A, 2)) :: temp1

        do count = 1, size(A, 1)
            test1 = A(count, :, :)
            call DGEEV('N', 'V', size(A, 2), test1, size(A, 2), temp1, dummy_val, dummy_vect, size(A, 2), temp2, &
                    size(A, 2), work, size(A, 2)*4, info)
            eig_val(count, :) = temp1
            eig_vect(count, :, :) = temp2
        end do

    end subroutine

    subroutine eig_complex(A, eig_val, eig_vect)
        ! Calculate eigenvalues and eigenvectors of 3D complex matrix A using LAPACK

        implicit none

        complex(kind=8), dimension(:, :, :), intent(in) :: A
        complex(kind=8), dimension(size(A, 1), size(A, 2)), intent(out) :: eig_val
        complex(kind=8), dimension(size(A, 1), size(A, 2), size(A, 3)), intent(out) :: eig_vect

        complex(kind=8) :: dummy(size(A, 2),size(A, 2)), work(size(A, 2)*2)
        complex(kind=8) :: Rwork(size(A, 2)*2)
        integer :: info, count

        complex(kind=8), dimension(size(A, 2), size(A, 2)) :: test1, temp2
        complex(kind=8), dimension(size(A, 2)) :: temp1

        do count = 1, size(A, 1)
            test1 = A(count, :, :)
            call ZGEEV('N', 'V', size(A, 2), test1, size(A, 2), temp1, dummy, size(A, 2), temp2, size(A, 2), work, &
                    size(A, 2)*2, Rwork, info)
            eig_val(count, :) = temp1
            eig_vect(count, :, :) = temp2
        end do

    end subroutine

end module functions