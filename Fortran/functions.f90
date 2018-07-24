module functions

    ! This module contains an assortment of common Fortran functions and subroutines.

contains

    function kron_rmat_eye(mat, val) result(B)
        ! Calculates np.kron(mat, np.eye(val)) of real square matrix mat
        ! Independent processes so number of OMP threads can take any value

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), intent(in) :: mat(:, :)
        integer, intent(in) :: val

        real(wp) :: temp(val, size(mat, 1), val, size(mat, 1))
        real(wp) :: B(size(mat, 1) * val, size(mat, 1) * val)
        integer count

        temp = 0
        do count = 1, val
            temp(count, :, count, :) = mat
        end do

        B = reshape(temp, [size(mat, 1) * val, size(mat, 1) * val])

    end function kron_rmat_eye

    function kron_cmat_eye(mat, val) result(B)
        ! Calculates np.kron(mat, np.eye(val)) of real square matrix mat
        ! Independent processes so number of OMP threads can take any value

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), intent(in) :: mat(:, :)
        integer, intent(in) :: val

        complex(wp) :: temp(val, size(mat, 1), val, size(mat, 1))
        complex(wp) :: B(size(mat, 1) * val, size(mat, 1) * val)
        integer count

        temp = 0
        do count = 1, val
            temp(count, :, count, :) = mat
        end do

        B = reshape(temp, [size(mat, 1) * val, size(mat, 1) * val])

    end function kron_cmat_eye

    function kron_eye_rmat(val, mat) result(B)
        ! Calculates np.kron(np.eye(val), mat) of real square matrix mat
        ! Independent processes so number of OMP threads can take any value

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), intent(in) :: mat(:, :)
        integer, intent(in) :: val

        real(wp) :: temp(size(mat, 1), val, size(mat, 1), val)
        real(wp) :: B(size(mat, 1) * val, size(mat, 1) * val)
        integer count

        temp = 0
        do count = 1, val
            temp(:, count, :, count) = mat
        end do

        B = reshape(temp, [size(mat, 1) * val, size(mat, 1) * val])

    end function kron_eye_rmat

    function kron_eye_cmat(val, mat) result(B)
        ! Calculates np.kron(np.eye(val), mat) of real square matrix mat
        ! Independent processes so number of OMP threads can take any value

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), intent(in) :: mat(:, :)
        integer, intent(in) :: val

        complex(wp) :: temp(size(mat, 1), val, size(mat, 1), val)
        complex(wp) :: B(size(mat, 1) * val, size(mat, 1) * val)
        integer count

        temp = 0
        do count = 1, val
            temp(:, count, :, count) = mat
        end do

        B = reshape(temp, [size(mat, 1) * val, size(mat, 1) * val])

    end function kron_eye_cmat

    function eye(A) result(B)
        ! Creates A x A identit matrix
        ! Independent processes so number of OMP threads can take any value

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        integer, intent(in) :: A
        real(wp) :: B(A, A)
        integer :: count

        B = 0._wp
        do count = 1, A
            B(count, count) = 1._wp
        end do

    end function eye

    function argsort(A) result(B)
        ! Calculates sorting indices of array A
        ! Iterative process so number of OMP threads must equal 1

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), dimension (:), intent(in) :: A
        real(wp) :: temp(size(A))
        integer(wp) :: index(1), B(size(A)), count

        temp =  A
        do count = 1, size(A)
            index = maxloc(temp)
            B(count) = index(1)
            temp(index) = -1E10_wp
        end do

    end function argsort

    function trace_real(A) result(C)
        ! Calculate trace of real matrix A
        ! Iterative process so number of OMP threads must equal 1

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), dimension (:, :), intent(in) :: A
        real(wp) :: C
        integer :: i
        C = 0

        do i = 1, size(A, 1)
            C = C + A(i, i)
        end do

    end function trace_real

    function trace_complex(A) result(C)
        ! Calculate trace of complex matrix A
        ! Iterative process so number of OMP threads must equal 1

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), dimension (:, :), intent(in) :: A
        complex(wp) :: C
        integer :: i
        C = 0

        do i = 1, size(A, 1)
            C = C + A(i, i)
        end do

    end function trace_complex

    function kron_real(A, B) result(AB)
        ! Calculate Kronecker product of real matrices A and B, adapted from Rosetta Code
        ! Iterative process so number of OMP threads must equal 1

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), intent(in) :: A(:, :), B(:, :)
        real(wp) :: AB(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2) )
        integer :: R, RA, RB, C, CA, CB, I, J

        R = 0.
        RA = size(A, DIM = 1)
        CA = size(A, DIM = 2)
        RB = size(B, DIM = 1)
        CB = size(B, DIM = 2)

        do I = 1, RA
            C = 0
            do J = 1, CA
                AB(R + 1 : R + RB, C + 1 : C + CB) = A(I, J) * B
                C = C + CB
            end do
            R = R + RB
        end do

    end function kron_real

    function kron_complex(A, B) result(AB)
        ! Calculate Kronecker product of complex matrices A and B, adapted from Rosetta Code
        ! Iterative process so number of OMP threads must equal 1

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), intent(in) :: A(:, :), B(:, :)
        complex(wp) :: AB(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2) )
        integer :: R, RA, RB, C, CA, CB, I, J

        R = 0
        RA = size(A, DIM = 1)
        CA = size(A, DIM = 2)
        RB = size(B, DIM = 1)
        CB = size(B, DIM = 2)

        do I = 1, RA
            C = 0
            do J = 1, CA
                AB(R + 1 : R + RB, C + 1 : C + CB) = A(I, J) * B
                C = C + CB
            end do
            R = R + RB
        end do

    end function kron_complex

    function inverse_real(A) result(B)
        ! Calculate inverse of 3D real matrix A using LAPACK
        ! Independent processes so number of OMP threads can take any value

        use omp_lib
        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), dimension(:, :, :), intent(in) :: A
        real(wp), dimension(size(A, 1), size(A, 2), size(A, 3)) :: B
        real(wp), dimension(size(A, 2), size(A, 3)) :: temp
        real(wp) :: work(size(A, 2)*2)

        integer :: count
        integer :: ipiv(size(A, 2)), info

        !!$omp parallel do default(private) &
        !!$omp& shared(A, B)
        do count = 1, size(A, 1)
            temp = A(count, :, :)
            call DGETRF(size(A, 2), size(A, 2), temp, size(A, 2), ipiv, info)
            call DGETRI(size(A, 2), temp, size(A, 2), ipiv, work, size(A, 2)*2, info)
            B(count, :, :) = temp
        end do
        !!$omp end parallel do

    end function inverse_real

    function inverse_complex(A) result(B)
        ! Calculate inverse of 3D complex matrix A using LAPACK
        ! Independent processes so number of OMP threads can take any value

        use omp_lib
        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), dimension(:, :, :), intent(in) :: A
        complex(wp), dimension(size(A, 1), size(A, 2), size(A, 3)) :: B
        complex(wp), dimension(size(A, 2), size(A, 3)) :: temp
        complex(wp) :: work(size(A, 2)*2)

        integer :: count
        integer :: ipiv(4), info

        !!$omp parallel do default(private) &
        !!$omp& shared(A, B)
        do count = 1, size(A, 1)
            temp = A(count, :, :)
            call ZGETRF(size(A, 2), size(A, 2), temp, size(A, 2), ipiv, info)
            call ZGETRI(size(A, 2), temp, size(A, 2), ipiv, work, size(A, 1)*2, info)
            B(count, :, :) = temp
        end do
        !!$omp end parallel do

    end function inverse_complex

    function expm_complex(A) result(B)
        ! Calculate matrix exponential of complex matrix A using Expokit

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), dimension(:, :), intent(in) :: A
        complex(wp), dimension(size(A, 1), size(A, 2)) :: B

        integer, parameter :: ideg = 2 ! Pade approximation, 6 is reccomended but 2 appears to be stable
        complex(wp) :: t = 1._wp
        complex(wp), dimension(4 * size(A, 1) * size(A, 2) + ideg + 1) :: wsp
        integer, dimension(size(A, 1)) :: iwsp
        integer :: iexp, ns, iflag

        call ZGPADM(ideg, size(A, 1), t, A, size(A, 1), wsp, size(wsp, 1), iwsp, iexp, ns, iflag)
        B = reshape(wsp(iexp : iexp + size(A, 1) * size(A, 1) - 1), [size(A, 1), size(A, 2)])

    end function expm_complex

    subroutine eig_real(A, eig_val, eig_vect)
        ! Calculate eigenvalues and eigenvectors of 3D real matrix A using LAPACK
        ! Independent processes so number of OMP threads can take any value

        use omp_lib
        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), dimension(:, :, :), intent(in) :: A
        real(wp), dimension(size(A, 1), size(A, 2)), intent(out) :: eig_val
        real(wp), dimension(size(A, 1), size(A, 2), size(A, 3)), intent(out) :: eig_vect

        real(wp) :: dummy_vect(size(A, 2), size(A, 2)), dummy_val(size(A, 2)), work(size(A, 2)*4)
        integer :: info, count

        real(wp), dimension(size(A, 2), size(A, 2)) :: test1, temp2
        real(wp), dimension(size(A, 2)) :: temp1

        !!$omp parallel do default(private) &
        !!$omp& shared(A, eig_val, eig_vect)
        do count = 1, size(A, 1)
            test1 = A(count, :, :)
            call DGEEV('N', 'V', size(A, 2), test1, size(A, 2), temp1, dummy_val, dummy_vect, size(A, 2), temp2, &
                    size(A, 2), work, size(A, 2)*4, info)
            eig_val(count, :) = temp1
            eig_vect(count, :, :) = temp2
        end do
        !!$omp end parallel do

    end subroutine

    subroutine eig_complex(A, eig_val, eig_vect)
        ! Calculate eigenvalues and eigenvectors of 3D complex matrix A using LAPACK.
        ! Independent processes so number of OMP threads can take any value

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        complex(wp), dimension(:, :, :), intent(in) :: A
        complex(wp), dimension(size(A, 1), size(A, 2)), intent(out) :: eig_val
        complex(wp), dimension(size(A, 1), size(A, 2), size(A, 3)), intent(out) :: eig_vect

        complex(wp) :: dummy(size(A, 2),size(A, 2)), work(size(A, 2)*2)
        complex(wp) :: Rwork(size(A, 2)*2)
        integer :: info, count

        complex(wp), dimension(size(A, 2), size(A, 2)) :: test1, temp2
        complex(wp), dimension(size(A, 2)) :: temp1

        !!$omp parallel do default(private) &
        !!$omp& shared(A, eig_val, eig_vect)
        do count = 1, size(A, 1)
            test1 = A(count, :, :)
            call ZGEEV('N', 'V', size(A, 2), test1, size(A, 2), temp1, dummy, size(A, 2), temp2, size(A, 2), work, &
                    size(A, 2)*2, Rwork, info)
            eig_val(count, :) = temp1
            eig_vect(count, :, :) = temp2
        end do
        !!$omp end parallel do

    end subroutine

end module functions