program example_main

    use example_subroutine
    implicit none

    real(kind = 8) :: eigval(100, 4), dummy(4, 4) !, work(8)

    real(kind = 8), dimension(4, 4) :: test1
    real(kind = 8), dimension(100, 4, 4) :: hamiltonian

    real(kind = 8) :: eig_vector_complex(100, 4, 4)
    integer :: info, count
    real(kind=8), dimension(4, 4) :: temp2
    real(kind=8), dimension(4) :: temp1

    real(kind=8) ::  dummy1(2), dummy2(4,4), work(4), work2(8), eig_vector(2, 2)
    double precision :: Rwork

    call calculate_hamiltonian(hamiltonian)
    write(6, *) 'hamiltonian(1, :, :)', hamiltonian(1, :, :)

    ! Calculate without temporary variables
!    do count = 1, size(hamiltonian, 1)
!        test1 = hamiltonian(count, :, :)
!!        call ZGEEV('N', 'V', 4, test1, 4, eigval(count, :), dummy, 4, &
!!                eig_vector_complex(count, :, :), 4, work, 8, Rwork, info)
!        call DGEEV('N', 'V', 2, test1, 2, eigval(count, :), &
!                dummy1, dummy2, 2, eig_vector_complex(count, :, :), 2, work2, 8, info)
!    end do

    ! Calculate with temporary variables
    do count = 1, size(hamiltonian, 1)
        test1 = hamiltonian(count, :, :)
!        call ZGEEV('N', 'V', 4, test1, 4, temp1, dummy, 4, &
!                temp2, 4, work, 8, Rwork, info)
        call DGEEV('N', 'V', 2, test1, 2, temp1, &
                dummy1, dummy2, 2, temp2, 2, work2, 8, info)
        eigval(count, :) = temp1
        eig_vector_complex(count, :, :) = temp2
    end do

    write(6, *) 'hamiltonian(1, :, :)', hamiltonian(1, :, :)
    write(6, *) 'eigval(1, :)', eigval(1, :)
    write(6, *) 'eig_vector_complex(1, :, :)', eig_vector_complex(1, :, :)

end program example_main