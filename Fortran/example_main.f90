program example_main

    use example_subroutine
    implicit none

    complex*16 :: eigval(2), dummy(2, 2), work(4), eig_vector(2, 2)
    double precision :: Rwork
    complex*16, dimension(2, 2) :: hamiltonian
    integer :: info

    eigval = 0
    dummy = 0
    work = 0
    eig_vector = 0
    Rwork = 0
    info = 0
    hamiltonian = 0

    call calculate_hamiltonian(hamiltonian)

    write(6, *) 'hamiltonian before', hamiltonian
    call ZGEEV('N', 'V', 2, hamiltonian, 2, eigval, dummy, 4, eig_vector, 2, work, 4, Rwork, info)

    write(6, *) 'hamiltonian after', hamiltonian
    write(6, *) 'eigval', eigval
    write(6, *) 'eig_vector', eig_vector
    write(6, *) 'info', info
    write(6, *) 'work', work

end program example_main
