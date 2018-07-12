program example_main

    use example_subroutine
    implicit none

    complex(kind = 8) :: eigval(2), dummy(2, 2), work(4), eig_vector(2, 2)
    complex(kind = 8) :: temp(2, 2), temp2(2, 2), temp1(2)
    real(kind = 8) :: Rwork
    complex(kind = 8), dimension(2, 2) :: hamiltonian
    integer :: info, count

    call calculate_hamiltonian(hamiltonian)

!    do count = 1, size(hamiltonian, 1)
!        temp = hamiltonian(count, :, :)
!        call ZGEEV('N', 'V', 2, temp, 2, eigval(count, :), dummy, 2, &
!                eig_vector(count, :, :), 2, work, 4, Rwork, info)
!    end do

    ! Calculate with temporary variables todo why does changing Rwork to complex give correct results? still breaks -O2
!    do count = 1, size(hamiltonian, 1)
!        temp = hamiltonian(count, :, :)
!        call ZGEEV('N', 'V', 2, temp, 2, temp1, dummy, 4, &
!                temp2, 2, work, 4, Rwork, info)
!        eigval(count, :) = temp1
!        eig_vector(count, :, :) = temp2
!    end do

    write(6, *) 'hamiltonian before', hamiltonian

    call ZGEEV('N', 'V', 2, hamiltonian, 2, eigval, dummy, 4, eig_vector, 2, work, 4, Rwork, info)

    ! Calculate with function
!    call eig_complex(hamiltonian, eigval, eig_vector)

!    write(6, *) 'hamiltonian(1, :, :)', hamiltonian(1, :, :)
!    write(6, *) 'eigval(1, :)', eigval(1, :)
!    write(6, *) 'eig_vector(1, :, :)', eig_vector(1, :, :)

    write(6, *) 'hamiltonian after', hamiltonian
    write(6, *) 'eigval', eigval
    write(6, *) 'eig_vector', eig_vector

end program example_main
