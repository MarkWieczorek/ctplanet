'''
Calculate the inertia tensor (to first order) given a radial density profile
and shape of all interfaces.
'''
import numpy as np
from scipy.linalg import eigh

import pyshtools as pyshtools


# ==== InertiaTensor ====

def InertiaTensor(hilm, rho, i_core, normalize=False, quiet=True):
    '''
    Calculate the inertia tensor given a radial density profile and shape of
    each interface.

    Usage
    -----
    I, A, B, C, M, R, angles = InertiaTensor(hilm, rho, i_core,
                                                [normalize, quiet])

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
    angles : ndarray, size(3,2)
        Matrix with each row containing the latitude and longitude
        coordinates (in degrees) of the principal moments A, B and C.

    Parameters
    ----------
    hilm : array of SHCoeffs class instances, size(i_core+1)
        Array of SHCoeffs class instances of the spherical harmonic
        coefficients of the relief at each interface. hilm[0] corresponds to
        r=0.
    rho : ndarray, size(i_core)
        Array of the densities of each layer, where index i corresponds to the
        density between interfaces i and i+1.
    i_core : int
        index corresponding to the top of the core.
    normalize : bool, optional, default = False
        If True, return all moments normalized by MR^2
    quiet : bool, optional, default = True
        If False, print additional information, including the locations of the
        axes of the principal moments and gravitational coefficients.
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

    if quiet is False:
        e = np.zeros((3, 2))
        e[0, 0] = 90. - np.rad2deg(np.arccos(vec[0, 2]))
        e[0, 1] = np.rad2deg(np.arctan2(vec[0, 1], vec[0, 0]))
        e[1, 0] = 90. - np.rad2deg(np.arccos(vec[1, 2]))
        e[1, 1] = np.rad2deg(np.arctan2(vec[1, 1], vec[1, 0]))
        e[2, 0] = 90. - np.rad2deg(np.arccos(vec[2, 2]))
        e[2, 1] = np.rad2deg(np.arctan2(vec[2, 1], vec[2, 0]))
        print('Mass core (kg) = {:e}'.format(mass))
        print('R core (m) = {:e}'.format(r_core))
        print('I / (MR^2) = ', II / mass / r_core**2)
        print('A / (MR^2) = {:e}'.format(A / mass / r_core**2))
        print('B / (MR^2) = {:e}'.format(B / mass / r_core**2))
        print('C / (MR^2) = {:e}'.format(C / mass / r_core**2))
        print('A (lat, lon) = ', e[0, 0], e[0, 1])
        print('B (lat, lon) = ', e[1, 0], e[1, 1])
        print('C (lat, lon) = ', e[2, 0], e[2, 1])
        print('a (m) = ', hilm[i_core].expand(lat=e[0, 0], lon=e[0, 1]))
        print('b (m) = ', hilm[i_core].expand(lat=e[1, 0], lon=e[1, 1]))
        print('c (m) = ', hilm[i_core].expand(lat=e[2, 0], lon=e[2, 1]))
        # print('\nC20 (unnorm) = ', -(II[2, 2] - (II[0, 0] + II[1, 1])/2.)
        #      / mass / r_core**2)
        # print('C21 (unnorm) = ', II[0, 2] / mass / r_core**2)
        # print('S21 (unnorm) = ', II[1, 2] / mass / r_core**2)
        # print('C22 (unnorm) = ', (II[1, 1] - II[0, 0]) / 4.
        #       / mass / r_core**2)
        # print('S22 (unnorm) = ', II[0, 1] / 2. / mass / r_core**2)

    if normalize:
        return II / mass / r_core**2, A / mass / r_core**2, \
            B / mass / r_core**2, C / mass / r_core**2, mass, r_core, vec
    else:
        return II, A, B, C, mass, r_core, vec
