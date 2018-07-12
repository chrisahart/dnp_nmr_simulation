program example_main

    use example_subroutine
    use iso_fortran_env
    implicit none

    integer, parameter :: size = 2 ** (2) ! Change (value) to change number of spins
    real(kind=real64) :: eigval(size), eig_vector(size, size), eig_vector_inv(size, size)
    real(kind=real64), dimension(size, size) :: hamiltonian

    ! Calculate Hamiltonian
    call calculate_hamiltonian(hamiltonian)

    ! Calculate eigenvalues and eigenvectors
    call eig(hamiltonian, eigval, eig_vector)

    ! Calculate inverse eigenvectors
    call inv(eig_vector, eig_vector_inv)

    ! Write eigenvalues and eigenvectors
    write(6, *) 'eigval', eigval
    write(6, *) 'eig_vector', eig_vector
    write(6, *) 'eig_vector_inv', eig_vector_inv

end program example_main
