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

    function kron_real(A, B) result(C)

        real(kind = 8), dimension (:, :), intent(in) :: A, B
        real(kind = 8), dimension (:, :), allocatable :: C
        integer :: i = 0, j = 0
        integer :: m = 0, n = 0, p = 0, q = 0
        allocate(C(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2)))

        C = 0

        write(6,*) (size(C, 1))
        write(6,*) (size(C, 2))
        write(6,*) (size(C))

        write(6,*) 'test', C(1:3, 1:3)

        do i = 1, size(A, 1)
            do j = 1, size(A, 2)

                write(6,*) 'loop', j

                n = (i - 1) * size(B, 1) + 1
                m = n + size(B, 1)
                p = (j - 1) * size(B, 2) + 1
                q = p + size(B, 2)

                write(6,*) 'n', n
                write(6,*) 'm', m
                write(6,*) 'p', p
                write(6,*) 'q', q

                C(1:3, 1:3) = A(i, j) * B
                write(6,*) 'once'
                C(n : m, p : q) = A(i, j) * B


            enddo
        enddo

    end function kron_real

    function kron_complex(A, B) result(C)

        complex(kind=8), dimension (:, :), intent(in) :: A, B
        complex(kind=8), dimension (:, :), allocatable :: C
        integer :: i = 0, j = 0
        integer :: m = 0, n = 0, p = 0, q = 0
        allocate(C(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2)))

        !!$omp parallel do default(private) &
        !!$omp& shared(A, B, C)
        do i = 1, size(A, 1)
            do j = 1, size(A, 2)
                n = (i - 1) * size(B, 1) + 1
                m = n + size(B, 1)
                p = (j - 1) * size(B, 2) + 1
                q = p + size(B, 2)
                C(n : m, p : q) = A(i, j) * B
            enddo
        enddo
        !!$omp end parallel do

    end function kron_complex

    function expm_real(t, H) result(expH)

        real(kind = 8), intent(in) :: t
        real(kind = 8), dimension(:, :), intent(in) :: H
        real(kind = 8), dimension(size(H, 1), size(H, 2)) :: expH

        integer, parameter :: ideg = 6
        real(kind = 8), dimension(4 * size(H, 1) * size(H, 2) + ideg + 1) :: wsp
        integer, dimension(size(H, 1)) :: iwsp
        integer :: iexp, ns, iflag, n

        n = size(H, 1)
        call DGPADM(ideg, n, t, H, n, wsp, size(wsp, 1), iwsp, iexp, ns, iflag)
        expH = reshape(wsp(iexp : iexp + n * n - 1), shape(expH))

        ! ideg, degre of the diagonal Pade to be used
        ! n the order of H
        ! t the timescale
        ! H the argument matrix
        ! wsp workspace variable
        ! iexp output
        ! ns number of scaling-squaring used
        ! iflag exit flag

    end function expm_real

    function expm_complex(t, H) result(expH)

        complex(kind=8), intent(in) :: t
        complex(kind=8), dimension(:, :), intent(in) :: H
        complex(kind=8), dimension(size(H, 1), size(H, 2)) :: expH

        integer, parameter :: ideg = 6
        complex(kind=8), dimension(4 * size(H, 1) * size(H, 2) + ideg + 1) :: wsp
        integer, dimension(size(H, 1)) :: iwsp
        integer :: iexp, ns, iflag, n

        n = size(H, 1)
        call ZGPADM(ideg, n, t, H, n, wsp, size(wsp, 1), iwsp, iexp, ns, iflag)
        expH = reshape(wsp(iexp : iexp + n * n - 1), shape(expH))

        ! ideg, degre of the diagonal Pade to be used
        ! n the order of H
        ! t the timescale
        ! H the argument matrix
        ! wsp workspace variable
        ! iexp output
        ! ns number of scaling-squaring used
        ! iflag exit flag

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

    subroutine eig_complex(hamil_complex, eigvals, eig_vector_complexs)

        !use omp_lib
        implicit none

        complex(kind=8), dimension(10000, 4, 4), intent(in) :: hamil_complex
        complex(kind=8), dimension(10000, 4), intent(out) :: eigvals
        complex(kind=8), dimension(10000, 4, 4), intent(out) :: eig_vector_complexs
        complex(kind=8) :: eig_vector_complex(10000, 4,4)

        complex(kind=8) :: dummy_eig(4,4), work_eig(8), temp_eig(4,4)
        real(kind=8) :: Rwork_eig
        integer :: info_eig, count

        write(6,*) 'start'

        do count = 1, size(hamil_complex, 1)

            temp_eig = hamil_complex(count, :, :)
!            write(6,*) temp_eig
!            write(6, *) eig_vector_complexs(count, 1, 1)
            call ZGEEV('N', 'V', 4, temp_eig, 4, eigvals(count, :), dummy_eig, 4,  &
                    eig_vector_complex(count, :, :), 4, work_eig, 8, Rwork_eig, info_eig)
!            write(6, *) eig_vector_complexs(count, 1, 1)
        end do

        write(6,*) 'finished'
        write(6, *) eig_vector_complex(1, :, :)
        write(6, *) eigvals(1, :)
        write(6, *) info_eig

        eig_vector_complexs = eig_vector_complex

    end subroutine

end module f2py_functions