Program LinearEquations2

! double precision kind constant
  integer, parameter :: dp = kind(1.d0)

  ! Calculate exp(t*H) for an N-by-N matrix H using Expokit.
  function expm(t, H) result(expH)
    real(dp), intent(in) :: t
    real(dp), dimension(:,:), intent(in) :: H
    real(dp), dimension(size(H,1),size(H,2)) :: expH

    ! Expokit variables
    external :: DGPADM
    integer, parameter :: ideg = 6
    real(dp), dimension(4*size(H,1)*size(H,2) + ideg + 1) :: wsp
    integer, dimension(size(H,1))  :: iwsp
    integer :: iexp, ns, iflag, n

    if (size(H,1) /= size(H,2)) then
       stop 'expm: matrix must be square'
    end if

    n = size(H,1)
    call DGPADM(ideg, n, t, H, n, wsp, size(wsp,1), iwsp, iexp, ns, iflag)
    expH = reshape(wsp(iexp:iexp+n*n-1), shape(expH))
  end function expm

end program LinearEquations2
