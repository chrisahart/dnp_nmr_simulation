module interactions

    ! This module contains subroutines for calculating interactions.

contains

    subroutine anisotropy_coefficients(electron_frequency, gtensor, orientation_se, c0, c1, c2, c3, c4)
        ! Calculate time independent g-anisotropy coefficients.

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), intent(in) :: electron_frequency, gtensor(3), orientation_se(3)
        real(wp) :: gx, gy, gz, ca, cb, cg, sa, sb, sg, r11, r12, r13, r21, r22, r23, r31, r32, r33
        real(wp), intent(out) :: c0, c1, c2, c3, c4

        gx = electron_frequency * gtensor(1)
        gy = electron_frequency * gtensor(2)
        gz = electron_frequency * gtensor(3)

        ca = cos(orientation_se(1))
        cb = cos(orientation_se(2))
        cg = cos(orientation_se(3))
        sa = sin(orientation_se(1))
        sb = sin(orientation_se(2))
        sg = sin(orientation_se(3))

        r11 = ca * cb * cg - sa * sg
        r12 = sa * cb * cg + ca * sg
        r13 = -sb * cg
        r21 = -ca * cb * sg - sa * cg
        r22 = -sa * cb * sg + ca * cg
        r23 = sb * sg
        r31 = ca * sb
        r32 = sa * sb
        r33 = cb

        c0 = (1._wp / 3._wp) * (gx + gy + gz)
        c1 = 2._wp * 2._wp ** (1._wp / 2._wp) / 3._wp * (gx * r11 * r31 + gy * r12 * r32 + gz * r13 * r33)
        c2 = 2._wp * 2._wp ** (1._wp / 2._wp) / 3._wp * (gx * r21 * r31 + gy * r22 * r32 + gz * r23 * r33)
        c3 = 1._wp / 3._wp * (gx * (r11 ** 2._wp - r21 ** 2._wp) + gy * (r12 ** 2._wp - r22 ** 2._wp) + &
                gz * (r13 ** 2 - r23 ** 2))
        c4 = 2._wp / 3._wp * (gx * r11 * r21 + gy * r22 * r12 + gz * r13 * r23)

    end subroutine

    subroutine anisotropy(c0, c1, c2, c3, c4, freq_rotor, electron_frequency, time, gtensor, ganisotropy)
        ! Calculate time dependent g-anisotropy.

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), parameter :: PI = 4._wp * ATAN(1._wp)
        real(wp), intent(in) :: freq_rotor, time, gtensor(3), electron_frequency
        real(wp), intent(in) :: c0, c1, c2, c3, c4
        real(wp), intent(out) :: ganisotropy
        real(wp) :: gx, gy, gz

        gx = electron_frequency * gtensor(1)
        gy = electron_frequency * gtensor(2)
        gz = electron_frequency * gtensor(3)

        ganisotropy = c0 + c1 * cos(2._wp * PI * freq_rotor * time) + &
                c2 * sin(2._wp * PI * freq_rotor * time) + &
                c3 * cos(4._wp * PI * freq_rotor * time) + &
                c4 * sin(4._wp * PI * freq_rotor * time)

    end subroutine anisotropy

    subroutine hyperfine(hyperfine_coupling, hyperfine_angles, freq_rotor, time, hyperfine_zz, hyperfine_zx)
        ! Calculate time dependent hyperfine.

        use iso_fortran_env
        implicit none

        integer, parameter :: wp = selected_real_kind(15, 307)
        real(wp), parameter :: PI = 4._wp * ATAN(1._wp)
        real(wp), intent(in) :: hyperfine_coupling, freq_rotor, time, hyperfine_angles(3)
        real(wp), intent(out) :: hyperfine_zz, hyperfine_zx

        hyperfine_zz = hyperfine_coupling * (-0.5_wp * (sin(hyperfine_angles(2)) ** 2._wp) * &
                cos(2._wp * (2._wp * PI * freq_rotor * time + &
                        hyperfine_angles(3))) + 2._wp ** 2._wp * sin(hyperfine_angles(2)) * &
                cos(hyperfine_angles(2)) * &
                cos(2 * PI * freq_rotor * time + hyperfine_angles(3)))

        hyperfine_zx = hyperfine_coupling * (-0.5_wp * sin(hyperfine_angles(2)) * cos(hyperfine_angles(2)) * &
                cos(2._wp * PI * freq_rotor * time + &
                        hyperfine_angles(3)) - (2._wp ** 0.5_wp / 4._wp) * (sin(hyperfine_angles(2)) ** 2._wp) * &
                cos(2._wp * (2._wp * PI * freq_rotor * time + hyperfine_angles(3))) + &
                (2._wp ** 0.5_wp / 4._wp) * (3._wp * (cos(hyperfine_angles(2)) ** 2._wp) - 1._wp))

    end subroutine hyperfine

end module interactions