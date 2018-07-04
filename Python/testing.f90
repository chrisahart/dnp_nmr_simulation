module testing

contains

    subroutine testing_2

        use iso_fortran_env
        use omp_lib
        implicit none

        integer, parameter :: WP = REAL64
        complex(kind = 8), parameter :: i = (0, 1)
        real(kind = 8), parameter :: PI = 4.D0 * DATAN(1.D0), Planck = 6.62607004E-34, Boltzmann = 1.38064852E-23

        integer :: count
        complex(kind = 8), dimension(4, 4) :: spin2_s_x, spin2_s_y, spin2_s_z, identity_size4, test
        complex(kind = 8), dimension(2, 2) :: spin_x, spin_y, spin_z, identity_size2

        ! Identity matrix
        identity_size2 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_size2)))


    end subroutine testing_2

!    function kron(A, B) result(C)
!
!        complex(kind = 8), dimension (:, :), intent(in) :: A, B
!        complex(kind = 8), dimension (:, :), allocatable :: C
!        integer :: i = 0, j = 0, k = 0, l = 0
!        integer :: m = 0, n = 0, p = 0, q = 0
!
!        allocate(C(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2)))
!        C = 0
!
!        do i = 1, size(A, 1)
!            do j = 1, size(A, 2)
!                n = (i - 1) * size(B, 1) + 1
!                m = n + size(B, 1)
!                p = (j - 1) * size(B, 2) + 1
!                q = p + size(B, 2)
!                C(n : m, p : q) = A(i, j) * B
!            enddo
!        enddo
!
!    end function kron
end module


!    test = expm(5, spin2_s_x)
!    write(*, *) test


!function expm(t, H) result(expH)
!
!    real, intent(in) :: t
!    real, dimension(:,:), intent(in) :: H
!    real, dimension(size(H,1),size(H,2)) :: expH
!
!    ! Expokit variables
!    external :: DGPADM
!    integer, parameter :: ideg = 6
!    real, dimension(4*size(H,1)*size(H,2) + ideg + 1) :: wsp
!    integer, dimension(size(H,1))  :: iwsp
!    integer :: iexp, ns, iflag, n
!
!    n = size(H,1)
!    call DGPADM(ideg, n, t, H, n, wsp, size(wsp,1), iwsp, iexp, ns, iflag)
!    expH = reshape(wsp(iexp:iexp+n*n-1), shape(expH))
!
!end function expm
