program testing

    use iso_fortran_env
    use omp_lib
    implicit none

    integer, parameter :: WP = REAL64
    complex(kind=8), parameter :: i = (0, 1)
    real(kind=8), parameter :: PI = 4.D0 * DATAN(1.D0), Planck = 6.62607004E-34, Boltzmann = 1.38064852E-23

    integer :: count
    complex(kind=8), dimension(4, 4) :: spin2_s_x, spin2_s_y, spin2_s_z, identity_size4, test
    complex(kind=8), dimension(2, 2) :: spin_x, spin_y, spin_z, identity_size2

    real :: t
    real, dimension(2,2) :: H
    real, dimension(size(H,1),size(H,2)) :: expH

    !external :: DGPADM
    integer, parameter :: ideg = 6
    real, dimension(4*size(H,1)*size(H,2) + ideg + 1) :: wsp
    integer, dimension(size(H,1))  :: iwsp
    integer :: iexp, ns, iflag, n

    ! Identity matrix
    identity_size2 = transpose(reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, 1.0_WP/), shape(identity_size2)))

    ! Pauli matrices
    spin_x = 0.5 * (reshape((/ 0.0_WP, 1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_x), order = (/2, 1/)))
    spin_y = 0.5 * i * (reshape((/ 0.0_WP, -1.0_WP, 1.0_WP, 0.0_WP/), shape(spin_y), order = (/2, 1/)))
    spin_z = 0.5 * (reshape((/ 1.0_WP, 0.0_WP, 0.0_WP, -1.0_WP/), shape(spin_z), order = (/2, 1/)))

    ! Expokit variables
    H = spin_x
    t = 5

    n = size(H,1)
    call DGPADM(ideg, n, t, H, n, wsp, size(wsp,1), iwsp, iexp, ns, iflag)
    expH = reshape(wsp(iexp:iexp+n*n-1), shape(expH))
    write(*, *) test

end program testing

function expm(t, H) result(expH)



end function expm
