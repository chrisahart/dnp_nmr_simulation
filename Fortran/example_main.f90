program example_main

    use example_subroutine
    use iso_fortran_env
    implicit none

!    real(kind=real64) :: eigval(2), dummy1(2), dummy2(4,4), work2(8), eig_vector(2, 2)
!    double precision :: Rwork
!    real(kind=real64), dimension(2, 2) :: hamiltonian
!    integer :: info
!
!    call calculate_hamiltonian(hamiltonian)
!
!    write(6, *) 'hamiltonian before', hamiltonian
!    call DGEEV('N', 'V', 2, hamiltonian, 2, eigval, dummy1, dummy2, 2, eig_vector, 2, work2, 8, info)
!
!    write(6, *) 'hamiltonian after', hamiltonian
!    write(6, *) 'eigval', eigval
!    write(6, *) 'eig_vector', eig_vector
!    write(6, *) 'info', info

    real(kind=real64) :: eigval(4), eig_vector(4, 4), eig_vector_inv(4, 4)
    real(kind=real64), dimension(4, 4) :: hamiltonian

    ! Calculate Hamiltonian
    call calculate_hamiltonian(hamiltonian)

    ! Calculate eigenvalues and eigenvectors
    call eig(hamiltonian, eigval, eig_vector)

    ! Calculate inverse eigenvectors
    call inv(eig_vector, eig_vector_inv)

    ! Write eigenvalues and eigenvectors
    write(6, *) 'hamiltonian', hamiltonian
    write(6, *) 'eigval', eigval
    write(6, *) 'eig_vector', eig_vector
    write(6, *) 'eig_vector_inv', eig_vector_inv

end program example_main
