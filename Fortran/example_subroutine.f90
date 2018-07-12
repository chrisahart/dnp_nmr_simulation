module example_subroutine

contains

    subroutine calculate_hamiltonian(hamiltonian)

        implicit none

        integer :: count
        complex(kind = 8), dimension(2, 2), intent(out) :: hamiltonian
        complex(kind = 8), dimension(2, 2) :: spin_x, spin_z

        spin_x = 0.5 * (reshape((/ 0.D0, 1.D0, 1.D0, 0.D0/), shape(spin_x), order = (/2, 1/)))
        spin_z = 0.5 * (reshape((/ 1.D0, 0.D0, 0.D0, -1.D0/), shape(spin_z), order = (/2, 1/)))

        hamiltonian =  2D6 * spin_z + 1D6 * spin_x + 1E6 * matmul(spin_x, spin_z)

!        do count = 1, 100
!            hamiltonian(count, :, :) =  2D6 * spin_z + 1D6 * spin_x + 1E6 * matmul(spin_x, spin_z)
!        end do

        write(6, *) 'hamiltonian', hamiltonian

    end subroutine calculate_hamiltonian

    subroutine eig_complex(hamiltonian_complex, eigval, eig_vector_complex)

        implicit none

        complex(kind=8), dimension(100, 2, 2), intent(in) :: hamiltonian_complex
        complex(kind=8), dimension(100, 2), intent(out) :: eigval
        complex(kind=8), dimension(100, 2, 2), intent(out) :: eig_vector_complex

        complex(kind=8) :: dummy(2,2), work(4)
        complex(kind=8) :: Rwork
        integer :: info, count

        complex(kind=8), dimension(2, 2) :: test1, temp2
        complex(kind=8), dimension(2) :: temp1

        ! Calculate energy, eigenvalues and eigenvectors of intrinsic Hamiltonian
        do count = 1, size(hamiltonian_complex, 1)
            test1 = hamiltonian_complex(count, :, :)
            call ZGEEV('N', 'V', 2, test1, 2, temp1, dummy, 2, temp2, 2, work, 4, Rwork, info)
            eigval(count, :) = temp1
            eig_vector_complex(count, :, :) = temp2
        end do

    end subroutine

end module