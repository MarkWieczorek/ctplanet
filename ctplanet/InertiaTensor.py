'''
Calculate the inertia tensor (to first order) given a radial density profile
and shape of all interfaces.
'''
import numpy as np
from scipy.linalg import eigh


# ==== InertiaTensor ====

def InertiaTensor_from_shape(hilm, rho, n, normalize=False, quiet=True):
    '''
    Calculate the inertia tensor given a radial density profile and shape of
    each interface.

    Returns
    -------
    I : ndarray, size(3, 3)
        The inertia tensor.
    A, B, C : float
        The principal moments of inertia, with A<B<C.
    M : float
        Mass of material beneath interface n.
    R : float
        The radius of the uppermost interface.
    angles : ndarray, size(3,2)
        Matrix with each row containing the latitude and longitude
        coordinates (in degrees) of the principal moments A, B and C.

    Parameters
    ----------
    hilm : array of SHCoeffs class instances, size(n+1)
        Array of SHCoeffs class instances of the spherical harmonic
        coefficients of the relief at each interface. hilm[0] corresponds to
        r=0 and hilm[n] to the uppermost interface.
    rho : ndarray, size(n)
        Array of the densities of each layer, where index i corresponds to the
        density between interfaces i and i+1.
    n : int
        index corresponding to the uppermost layer.
    normalize : bool, optional, default = False
        If True, return all moments normalized by MR^2
    quiet : bool, optional, default = True
        If False, print additional information, including the directions of the
        axes of the principal moments and gravitational coefficients.
    '''

    # Determine the contribution to the mass, average moment of inertia
    # and components of the inertia tensor for each layer.
    I0 = 0
    mass = 0
    r_n = hilm[n].coeffs[0, 0, 0]
    II = np.zeros((3, 3))

    for i in range(0, n):
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
        print('Mass (kg) = {:e}'.format(mass))
        print('R (m) = {:e}'.format(r_n))
        print('I / (MR^2) = ', II / mass / r_n**2)
        print('A / (MR^2) = {:e}'.format(A / mass / r_n**2))
        print('B / (MR^2) = {:e}'.format(B / mass / r_n**2))
        print('C / (MR^2) = {:e}'.format(C / mass / r_n**2))
        print('A (lat, lon) = ', e[0, 0], e[0, 1])
        print('B (lat, lon) = ', e[1, 0], e[1, 1])
        print('C (lat, lon) = ', e[2, 0], e[2, 1])
        print('a (m) = ', hilm[n].expand(lat=e[0, 0], lon=e[0, 1]))
        print('b (m) = ', hilm[n].expand(lat=e[1, 0], lon=e[1, 1]))
        print('c (m) = ', hilm[n].expand(lat=e[2, 0], lon=e[2, 1]))
        # print('\nC20 (unnorm) = ', -(II[2, 2] - (II[0, 0] + II[1, 1])/2.)
        #      / mass / r_n**2)
        # print('C21 (unnorm) = ', II[0, 2] / mass / r_n**2)
        # print('S21 (unnorm) = ', II[1, 2] / mass / r_n**2)
        # print('C22 (unnorm) = ', (II[1, 1] - II[0, 0]) / 4.
        #       / mass / r_n**2)
        # print('S22 (unnorm) = ', II[0, 1] / 2. / mass / r_n**2)

    if normalize:
        return II / mass / r_n**2, A / mass / r_n**2, \
            B / mass / r_n**2, C / mass / r_n**2, mass, r_n, vec
    else:
        return II, A, B, C, mass, r_n, vec


# === Compute the three moments of inertia given C and the gravity coefficients

def InertiaTensor_from_C(C, potential, normalize=False, r_norm=None,
                         quiet=True):
    '''
    Calculate the inertia tensor given the polar moment of inertia and the
    gravitational potential coefficients.

    Returns
    -------
    I : ndarray, size(3, 3)
        The inertia tensor.
    A, B, C : float
        The principal moments of inertia, with A<B<C.
    angles : ndarray, size(3,2)
        Matrix with each row containing the latitude and longitude
        coordinates (in degrees) of the principal moments A, B and C.

    Parameters
    ----------
    C : float
        The polar moment of inertia, which is assumed to be equal to the
        I33 component of the inertia tensor.
    potential : SHGravCoeffs
        An SHGravCoeffs instance containing the gravitational potential
        coefficients.
    normalize : bool, optional, default = False
        If True, return all moments normalized by MR^2
    r_norm : float, optional, default = None
        If specified, and if normalize is True, use this radius to normalize
        all output moments of inertia. If normalize is False, then r_norm will
        be used when printing the normalized moments to screen when quiet is
        False.
    quiet : bool, optional, default = True
        If False, print additional information, including the principal and
        normalized principal moments of inertial, and the directions of the
        axes of the principal moments.

    Notes
    -----
    This routine assumes that the polar moment of inertia C is equal to the
    I33 term of the inertia tensor. This is equivalent to assuming that the
    coordinate system defining the gravitational potential is aligned with
    the principal moment C. As such, the gravitational potential terms of order
    2 and degree 1 should be identically zero. If they are not, the returned
    value of the largest principal moment will differ slightly from the input
    value, and the difference provides an estimate of the error associated with
    the assumption that C=I33.
    '''

    if r_norm is None:
        r_norm = potential.r0

    mass = potential.mass
    r0 = potential.r0

    clm_unnorm = potential.to_array(normalization='unnorm', csphase=1, lmax=2,
                                    errors=False)

    I33 = C
    I22 = mass * r0**2 * (clm_unnorm[0, 2, 0] + 2 * clm_unnorm[0, 2, 2]) \
        + I33
    I11 = mass * r0**2 * (clm_unnorm[0, 2, 0] - 2 * clm_unnorm[0, 2, 2]) \
        + I33
    I12 = - 2 * clm_unnorm[1, 2, 2] * mass * r0**2
    I13 = - clm_unnorm[0, 2, 1] * mass * r0**2
    I23 = - clm_unnorm[1, 2, 1] * mass * r0**2

    II = np.array([[I11, I12, I13], [I12, I22, I23], [I13, I23, I33]])

    eig, vec = eigh(II)
    AA = eig[0]
    BB = eig[1]
    CC = eig[2]

    if quiet is False:
        e = np.zeros((3, 2))
        e[0, 0] = 90. - np.rad2deg(np.arccos(vec[0, 2]))
        e[0, 1] = np.rad2deg(np.arctan2(vec[0, 1], vec[0, 0]))
        e[1, 0] = 90. - np.rad2deg(np.arccos(vec[1, 2]))
        e[1, 1] = np.rad2deg(np.arctan2(vec[1, 1], vec[1, 0]))
        e[2, 0] = 90. - np.rad2deg(np.arccos(vec[2, 2]))
        e[2, 1] = np.rad2deg(np.arctan2(vec[2, 1], vec[2, 0]))

        print('I / (MR^2) = ', II / mass / r_norm**2)
        print('\nA / (MR^2) = {:e}'.format(AA / mass / r_norm**2))
        print('B / (MR^2) = {:e}'.format(BB / mass / r_norm**2))
        print('C / (MR^2) = {:e}'.format(CC / mass / r_norm**2))
        print('I_ave / (MR^2) = {:e}\n'.format(II.trace() / 3 / mass /
                                               r_norm**2))
        print('A  = {:e}'.format(AA))
        print('B = {:e}'.format(BB))
        print('C = {:e}'.format(CC))
        print('I_ave = {:e}\n'.format(II.trace() / 3))
        print('A (lat, lon) = ', e[0, 0], e[0, 1])
        print('B (lat, lon) = ', e[1, 0], e[1, 1])
        print('C (lat, lon) = ', e[2, 0], e[2, 1])

    if normalize:
        return II / mass / r_norm**2, AA / mass / r_norm**2, \
            BB / mass / r_norm**2, CC / mass / r_norm**2, vec
    else:
        return II, AA, BB, CC, vec


# === moi : calculate the mean normalized moment of inertia

def moi(radius, rho, n, normalized=True):
    """
    Calculate the mean, normalized, moment of inertia up to index n.

    Returns
    -------
    x : float
        The mean moment of inertia computed using a 1-D density profile. If
        normalized is True, the moment of inertial will be normalized by
        M R^2.

    Parameters
    ----------
    radius : ndarray, size (n+1)
        A vector of radii, where radius[0] is the center of the planet and
        radius[n] is the surface.
    rho : ndarray, size (n+1)
        A vector of densities of the layers, where rho[i] is the density
        between radius[i] and radius[i+1]. The density above the surface,
        rho[n], is set to zero. The density of the base of the crust is
        rho[i_crust], the density of the upper mantle is rho[i_crust-1], and
        the density if the upper core is rho[i_core-1].
    n : integer
        Maximum indice of the radius to use when computing the moment of
        inertia.
    normalized : bool, optional, default = True
        If True, return the mean moment of inertia normalized by M R^2, where
        R is radius[n] and M is the mass below radius[n].
    """
    II = 0.
    mass = 4. * np.pi / 3. * rho[0] * radius[1]**3

    for i in range(2, n+1):
        mass += 4. * np.pi / 3. * rho[i-1] * (radius[i]**3 - radius[i-1]**3)

        II += 8. * np.pi / 15. * rho[i-1] * (radius[i]**5 - radius[i-1]**5)

    if normalized:
        return II / mass / radius[n]**2
    else:
        return II
