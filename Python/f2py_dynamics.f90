module f2py_dynamics

contains

  subroutine spin2_electronuclear()

    use iso_fortran_env
    use Kronecker

    implicit none

    integer, parameter :: WP = REAL64
    complex*8, dimension(4, 4) :: spin2_s_x
    complex*8, dimension(2, 2) :: spin_s_x, spin_s_y, spin_s_z

    ! Pauli matrices
    spin_s_x = 0.5 * (reshape((/ 0.0_WP, 1.0_WP, &
                                  1.0_WP, 0.0_WP/), &
                                  shape(spin_s_x), order=(/2, 1/)))

    !call kronecker(spin_s_x, spin2_s_x)
    KronProd(spin_s_x, spin_s_x, spin2_s_x)
    write(6,*) spin2_s_x


    end subroutine spin2_electronuclear

    subroutine kronecker(spin_s_x, spin2_s_x)

      use iso_fortran_env

      implicit none

      integer, parameter :: WP = REAL64

      complex*8, dimension(2, 2), intent(in) :: spin_s_x
      complex*8, dimension(4, 4), intent(out) :: spin2_s_x

      ! 4x4 Matrix operators for S operator
      spin2_s_x = 0.5 * (reshape((/ 0.0_WP, 0.0_WP, 1.0_WP, 0.0_WP, &
                                    0.0_WP, 0.0_WP, 0.0_WP, 1.0_WP, &
                                    1.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
                                    0.0_WP, 1.0_WP,  0.0_WP, 0.0_WP/), &
                                    shape(spin2_s_x), order=(/2, 1/)))

    return
    end subroutine kronecker

  end module
