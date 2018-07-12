module interactions

contains

!    subroutine anisotropy_coefficients(electron_frequency, gtensor, orientation_se, c0, c1, c2, c3, c4)
!
!        implicit none
!
!        real(kind = 8), intent(in) :: electron_frequency, gtensor(3), orientation_se(3)
!        real(kind = 8) :: gx, gy, gz, ca, cb, cg, sa, sb, sg, r11, r12, r13, r21, r22, r23, r31, r32, r33
!        real(kind = 8), intent(out) :: c0, c1, c2, c3, c4
!
!        gx = electron_frequency * gtensor(1)
!        gy = electron_frequency * gtensor(2)
!        gz = electron_frequency * gtensor(3)
!
!        ca = cos(orientation_se(1))
!        cb = cos(orientation_se(2))
!        cg = cos(orientation_se(3))
!        sa = sin(orientation_se(1))
!        sb = sin(orientation_se(2))
!        sg = sin(orientation_se(3))
!
!        r11 = ca * cb * cg - sa * sg
!        r12 = sa * cb * cg + ca * sg
!        r13 = -sb * cg
!        r21 = -ca * cb * sg - sa * cg
!        r22 = -sa * cb * sg + ca * cg
!        r23 = sb * sg
!        r31 = ca * sb
!        r32 = sa * sb
!        r33 = cb
!
!        c0 = 1.D0 / 3.D0 * (gx + gy + gz)
!        c1 = 2.D0 * 2.D0 ** (1.D0 / 2.D0) / 3.D0 * (gx * r11 * r31 + gy * r12 * r32 + gz * r13 * r33)
!        c2 = 2.D0 ** (1.D0 / 2.D0) / 3.D0 * (gx * r21 * r31 + gy * r22 * r32 + gz * r23 * r33)
!        c3 = 1.D0 / 3.D0 * (gx * (r11 ** 2.D0 - r21 ** 2.D0) + gy * (r12 ** 2.D0 - r22 ** 2.D0) + &
!                gz * (r13 ** 2 - r23 ** 2))
!        c4 = 2.D0 / 3.D0 * (gx * r11 * r21 + gy * r22 * r12 + gz * r13 * r23)
!
!    end subroutine
!
!    subroutine anisotropy(c0, c1, c2, c3, c4, freq_rotor, time, g_anisotropy)
!
!        implicit none
!
!        real(kind = 8), parameter :: PI = 4.D0 * DATAN(1.D0)
!        real(kind = 8), intent(in) :: freq_rotor, time
!        real(kind = 8), intent(in) :: c0, c1, c2, c3, c4
!        real(kind = 8), intent(out) :: g_anisotropy
!
!        gx = electron_frequency * gtensor(1)
!        gy = electron_frequency * gtensor(2)
!        gz = electron_frequency * gtensor(3)
!
!        g_anisotropy = c0 + c1 * cos(2.D0 * PI * freq_rotor * time) + &
!                c2 * sin(2.D0 * PI * freq_rotor * time) + &
!                c3 * cos(4.D0 * PI * freq_rotor * time) + &
!                c4 * sin(4.D0 * PI * freq_rotor * time)
!
!    end subroutine anisotropy
!
!    subroutine hyperfine(hyperfine_coupling, hyperfine_angles, time, hyperfine_zz, hyperfine_zx)
!
!        implicit none
!
!        real(kind = 8), parameter :: PI = 4.D0 * DATAN(1.D0)
!        real(kind = 8), intent(in) :: hyperfine_coupling, time, hyperfine_angles(3)
!        real(kind = 8), intent(out) :: hyperfine_zz, hyperfine_zx
!
!        hyperfine_zz = hyperfine_coupling * (-0.5D0 * (sin(hyperfine_angles(2)) ** 2.D0) * &
!                cos(2.D0 * (2.D0 * PI * freq_rotor * time + &
!                        hyperfine_angles(3))) + 2.D0 ** 2.D0 * sin(hyperfine_angles(2)) * &
!                cos(hyperfine_angles(2)) * &
!                cos(2 * PI * freq_rotor * time + hyperfine_angles(3)))
!
!        hyperfine_zx = hyperfine_coupling * (-0.5D0 * sin(hyperfine_angles(2)) * cos(hyperfine_angles(2)) * &
!                cos(2.D0 * PI * freq_rotor * time + &
!                        hyperfine_angles(3)) - (2.D0 ** 0.5D0 / 4.D0) * (sin(hyperfine_angles(2)) ** 2.D0) * &
!                cos(2.D0 * (2.D0 * PI * freq_rotor * time + hyperfine_angles(3))) + &
!                (2.D0 ** 0.5D0 / 4.D0) * (3.D0 * (cos(hyperfine_angles(2)) ** 2.D0) - 1.D0))
!
!    end subroutine hyperfine

end module interactions