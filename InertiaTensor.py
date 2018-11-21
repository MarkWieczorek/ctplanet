'''
Calculate the inertia tensor (to first order) given a radial density profile
and shape of all interfaces.
'''
import numpy as np
from scipy.linalg import eigh

import pyshtools as pyshtools


# ==== InertiaTensor ====

def InertiaTensor(hilm, rho, i_core, normalize=False):
    '''
    Calculate the Inertia tensor given a radial density profile and shape of
    each interface.

    Usage
    -----
    I, A, B, C, M, R, I_vectors = InertiaTensor(hilm, rho, i_core)

    Returns
    -------
    I : ndarray, size(3, 3)
        The inertia tensor.
    A, B, C : float
        The principal moments of inertia, with A<B<C.
    M : float
        Mass of the core.
    R : float
        The core radius.
    I_vectors : ndarray, size(3,3)
        Matrix with each column containing the directions of the principal
        moments A, B and C.

    Parameters
    ----------
    hilm : array of SHCoeffs class instances, size(i_core+1)
        Array of SHCoeffs class instances of the spherical harmonic
        coefficients of the hydrostatic relief at each interface. hilm[0]
        corresponds to r=0.
    rho : ndarray, size(i_core)
        Array of the densities of each layer, where index i corresponds to the
        density between interfaces i and i+1.
    i_core : int
        index corresponding to the top of the core.
    normalize : bool, optional
        If True, return moments normalized by MR^2
    '''

    # Determine the contribution to the mass, average moment of inertia
    # and components of the inertia tensor for each layer.
    I0 = 0
    mass = 0
    r_core = hilm[i_core].coeffs[0, 0, 0]
    II = np.zeros((3, 3))

    for i in range(0, i_core):
        r1 = hilm[i].coeffs[0, 0, 0]
        r2 = hilm[i+1].coeffs[0, 0, 0]
        I0 += rho[i] * np.pi * 8. / 15. * (r2**5 - r1**5)
        mass += rho[i] * 4. * np.pi / 3. * (r2**3 - r1**3)

        # Ixx
        II[0, 0] += rho[i] * 4. * np.pi * (1. / (3. * np.sqrt(5)) *
                                           (r2**4 * hilm[i+1].coeffs[0, 2, 0] -
                                           r1**4 * hilm[i].coeffs[0, 2, 0])
                                           - np.sqrt(12. / 5.) / 6. *
                                           (r2**4 * hilm[i+1].coeffs[0, 2, 2] -
                                            r1**4 * hilm[i].coeffs[0, 2, 2])
                                           )

        # Iyy
        II[1, 1] += rho[i] * 4. * np.pi * (1. / (3. * np.sqrt(5)) *
                                           (r2**4 * hilm[i+1].coeffs[0, 2, 0] -
                                           r1**4 * hilm[i].coeffs[0, 2, 0])
                                           + np.sqrt(12. / 5.) / 6. *
                                           (r2**4 * hilm[i+1].coeffs[0, 2, 2] -
                                            r1**4 * hilm[i].coeffs[0, 2, 2])
                                           )

        # Izz
        II[2, 2] += rho[i] * 4. * np.pi * (2. / (3. * np.sqrt(5)) *
                                           (r1**4 * hilm[i].coeffs[0, 2, 0] -
                                           r2**4 * hilm[i+1].coeffs[0, 2, 0])
                                           )

        # Ixy
        II[1, 0] += rho[i] * 4. * np.pi * (np.sqrt(12. / 5.) / 6. *
                                           (r2**4 * hilm[i+1].coeffs[1, 2, 2] -
                                           r1**4 * hilm[i].coeffs[1, 2, 2])
                                           )

        # Iyz
        II[2, 1] += rho[i] * 4. * np.pi * (np.sqrt(3. / 5.) / 3. *
                                           (r2**4 * hilm[i+1].coeffs[1, 2, 1] -
                                           r1**4 * hilm[i].coeffs[1, 2, 1])
                                           )

        # Ixz
        II[2, 0] += rho[i] * 4. * np.pi * (np.sqrt(3. / 5.) / 3. *
                                           (r2**4 * hilm[i+1].coeffs[0, 2, 1] -
                                           r1**4 * hilm[i].coeffs[0, 2, 1])
                                           )

    II[0, 0] += I0
    II[1, 1] += I0
    II[2, 2] += I0
    II[0, 1] = II[1, 0]
    II[0, 2] = II[2, 0]
    II[1, 2] = II[2, 1]

    eig, vec = eigh(II)
    A = eig[0]
    B = eig[1]
    C = eig[2]

    print('A (lat, lon) = ', 90 - np.rad2deg(np.arccos(vec[0, 2])),
          np.rad2deg(np.arctan2(vec[0, 1], vec[0, 0])))
    print('B (lat, lon) = ', 90 - np.rad2deg(np.arccos(vec[1, 2])),
          np.rad2deg(np.arctan2(vec[1, 1], vec[1, 0])))
    print('C (lat, lon) = ', 90 - np.rad2deg(np.arccos(vec[2, 2])),
          np.rad2deg(np.arctan2(vec[2, 1], vec[2, 0])))

    print('C20 (unnorm) = ', -(C - (A + B)/2) / mass / r_core**2)
    print('C21 (unnorm) = ', II[0, 2] / mass / r_core**2)
    print('S21 (unnorm) = ', II[1, 2] / mass / r_core**2)
    print('C22 (unnorm) = ', (B - A) / 4 / mass / r_core**2)
    print('S22 (unnorm) = ', II[0, 1] / 2 / mass / r_core**2)

    if normalize:
        return II / mass / r_core**2, A / mass / r_core**2, \
            B / mass / r_core**2, C / mass / r_core**2, mass, r_core, vec
    else:
        return II, A, B, C, mass, r_core, vec
