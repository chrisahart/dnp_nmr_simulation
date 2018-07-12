module example_subroutine

contains

    subroutine calculate_hamiltonian(hamiltonian)

        implicit none

        complex*16, dimension(2, 2), intent(out) :: hamiltonian
        complex*16, dimension(2, 2) :: spin_x, spin_z

        hamiltonian = 0
        spin_x = 0
        spin_z = 0

!        spin_x = 0.5 * (reshape((/ 0.D0, 1.D0, 1.D0, 0.D0/), shape(spin_x), order = (/2, 1/)))
!        spin_z = 0.5 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))

        spin_x = 0.5 * (reshape((/ 0.D0, 1.D0, 1.D0, 0.D0/), shape(spin_x), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))

        hamiltonian =  2D6 * spin_z + 1D6 * spin_x + 1E6 * matmul(spin_x, spin_z)

        write(6, *) 'hamiltonian', hamiltonian

    end subroutine calculate_hamiltonian

end module