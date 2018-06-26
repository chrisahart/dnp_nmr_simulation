import numpy as np
import parameters as param


def anisotropy_coefficients(angles):
    """ Calculate time independent coefficients for electron g-anisotropy.
    """

    gx = param.electron_frequency * (2.00614 / 2) + 18.76e6
    gy = param.electron_frequency * (2.00194 / 2) + 92.4e6
    gz = param.electron_frequency * (2.00988 / 2) + 18.2e6

    ca = np.cos(angles[0])
    cb = np.cos(angles[1])
    cg = np.cos(angles[2])
    sa = np.sin(angles[0])
    sb = np.sin(angles[1])
    sg = np.sin(angles[2])

    r11 = ca * cb * cg - sa * sg
    r12 = sa * cb * cg + ca * sg
    r13 = -sb * cg
    r21 = -ca * cb * sg - sa * cg
    r22 = -sa * cb * sg + ca * cg
    r23 = sb * sg
    r31 = ca * sb
    r32 = sa * sb
    r33 = cb

    c0 = 1 / 3 * (gx + gy + gz)
    c1 = 2 * np.sqrt(2) / 3 * (gx * r11 * r31 + gy * r12 * r32 + gz * r13 * r33)
    c2 = 2 * np.sqrt(2) / 3 * (gx * r21 * r31 + gy * r22 * r32 + gz * r23 * r33)
    c3 = 1 / 3 * (gx * (r11 ** 2 - r21 ** 2) + gy * (r12 ** 2 - r22 ** 2) + gz * (r13 ** 2 - r23 ** 2))
    c4 = 2 / 3 * (gx * r11 * r21 + gy * r22 * r12 + gz * r13 * r23)

    return c0, c1, c2, c3, c4
