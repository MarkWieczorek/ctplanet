'''
Functions for calculating the shape of hydrostatic density interfaces and their
gravitational potential in a planet with a non-hydrostatic lithosphere.
'''
import numpy as np
import scipy.linalg.lapack as lapack

import pyshtools as pysh


# ==== HydrostaticShapeLith ====

def HydrostaticShapeLith(radius, rho, ilith, potential, topo, rho_surface,
                         r_sigma, omega, lmax, rp=None, mp=None, nmax=7):
    """
    Calculate the shape of hydrostatic relief in a rotating planet or moon with
    a non-hydrostatic lithosphere, along with the total gravitation potential
    of the hydrostatic interfaces. For the case of a moon in synchronous
    rotation, optionally include the tidal potential.

    Returns
    -------
    hlm : array of SHCoeffs class instances, size(n+1)
        Array of SHCoeffs class instances of the spherical harmonic
        coefficients of the hydrostatic relief at each interface.
    clm_hydro : SHCoeffs class instance containing the gravitational potential
        resulting from all hydrostatic interfaces.
    mass : float
        Total mass of the planet, assuming a spherical shape and the provided
        1D density profile.

    Parameters
    ----------
    radius : ndarray, float, size(n+1)
        Radius of each density interface, where index 0 corresponds to the
        center of the planet and n corresponds to the surface.
    density : ndarray, float, size(n+1)
        Density of layer i between radius[i] and radius[i+1]. The density
        at index 0 is from the center of the planet to radius[1], whereas the
        the density at index n should be zero.
    ilith : int
        Index of the interface that corresponds to the base of the lithosphere.
    potential : SHGravCoeffs class instance
        Observed gravitational potential coefficients.
    topo : SHCoeffs class instance
        Observed shape of the planet.
    rho_surface : float
        Effective density of the surface relief.
    r_sigma : float
        Radius of the mass sheet source in the lithosphere.
    omega : float
        Angular rotation rate of the planet.
    lmax : int
        Maximum spherical harmonic degree to compute for the hydrostatic
        relief at each interface.
    rp : float, optional, default = None
        If specified, include the tidal potential acting on a synchronously
        rotating moon, where rp is the average distance between the planet
        and satellite.
    mp : float, optional, default = None
        The mass of the host planet, at a distance rp from the satellite.
    nmax : int, optional, default = 7
        The order of the approximation when computing the gravitational
        potential of the surface.
    """
    tides = False
    if rp is not None:
        if mp is None:
            raise ValueError('When including tides, both rp and mp must be ' +
                             'specified.')
        tides = True
    if mp is not None:
        if rp is None:
            raise ValueError('When including tides, both rp and mp must be ' +
                             'specified.')
        tides = True

    if len(radius) != len(rho):
        raise ValueError('Length of radius and density must be the same.' +
                         'len(radius) = {:d}. len(density) = {:d}.'
                         .format(len(radius), len(rho)))

    n = len(radius) - 1  # index of surface
    lmaxgrid = 3*lmax  # increase grid size to avoid aliasing

    g = pysh.constants.G.value
    gm = potential.gm
    r_ref = potential.r0

    hlm = [pysh.SHCoeffs.from_zeros(lmax) for i in range(ilith+1)]
    clm_hydro = pysh.SHCoeffs.from_zeros(lmax)

    for i in range(ilith+1):
        hlm[i].coeffs[0, 0, 0] = radius[i]

    # First determine the spherical harmonic coefficients of (Y20 Ylm)
    # and for tides, (Y22 Ylm). We are only concerned with the coefficient
    # that corresponds to lm, so for each lm, store only the lm component in
    # the array cp20 and cp22.

    sh20 = np.zeros((2, lmax+1, lmax+1))
    sh22 = np.zeros((2, lmax+1, lmax+1))
    sh = np.zeros((2, lmax+1, lmax+1))
    cp20 = np.zeros((2, lmax+1, lmax+1))
    cp22 = np.zeros((2, lmax+1, lmax+1))

    sh20[0, 2, 0] = 1.  # Y20
    sh22[0, 2, 2] = 1.  # Y22

    for l in range(1, lmax+1):
        for m in range(0, l+1):
            sh[0, l, m] = 1.
            coeffs = pysh.expand.SHMultiply(sh20, sh)
            cp20[0, l, m] = coeffs[0, l, m]
            if m != 0:
                cp20[1, l, m] = cp20[0, l, m]
            if l == 2 and m == 0:
                p402020 = coeffs[0, 4, 0]
            if l == 2 and m == 1:
                p412021 = coeffs[0, 4, 1]
            if l == 2 and m == 2:
                p422022 = coeffs[0, 4, 2]

            coeffs = pysh.expand.SHMultiply(sh22, sh)
            cp22[0, l, m] = coeffs[0, l, m]
            sh[0, l, m] = 0.

            if m > 0:
                sh[1, l, m] = 1.
                coeffs = pysh.expand.SHMultiply(sh22, sh)
                cp22[1, l, m] = coeffs[1, l, m]
                sh[1, l, m] = 0.

    # Calculate delta_rho

    drho = np.zeros(n+1)
    mass = np.zeros(n+1)

    for i in range(1, n):
        drho[i] = rho[i-1] - rho[i]
    drho[n] = rho_surface  # Effective density of the surface layer

    a = np.zeros((ilith+2, ilith+2))
    atides = np.zeros((2, ilith+1))
    b4 = np.zeros((2, 3, ilith+1))

    # Calculate cumulate mass function

    for i in range(1, n+1):
        if i == 1:
            mass[1] = 4. * np.pi * radius[1]**3 * rho[0] / 3.
        else:
            mass[i] = mass[i-1] + 4. * np.pi * \
                (radius[i]**3 - radius[i-1]**3) * rho[i-1] / 3.

    mass_model = mass[n]

    # Calculate potential coefficients of the surface relief

    grid = topo.expand(grid='DH2', lmax=lmaxgrid, lmax_calc=lmax, extend=False)
    cplus, r_surface = pysh.gravmag.CilmPlusDH(grid.data, nmax, gm/g,
                                               rho_surface, lmax=lmax)
    cminus, r_surface = pysh.gravmag.CilmMinusDH(grid.data, nmax, gm/g,
                                                 rho_surface, lmax=lmax)

    # Calculate matrix A and invert for relief.

    for l in range(1, lmax+1):
        for m in range(0, lmax+1):
            for i in range(1, ilith+1):  # zero index not computed
                for j in range(1, ilith+1):
                    if i == j:  # cp20 for sin and cosine terms are equal
                        a[i, j] = 4. * np.pi * g * drho[i] * radius[i] / \
                            (2. * l + 1.) - g * mass[i] / radius[i]**2 + \
                            (2./3.) * radius[i] * omega**2 * \
                            (1. - cp20[0, l, m] / np.sqrt(5.0))
                    elif j < i:
                        a[i, j] = 4. * np.pi * g * drho[j] / (2. * l + 1.) \
                            * radius[j] * (radius[j] / radius[i])**(l+1)
                    else:
                        a[i, j] = 4. * np.pi * g * drho[j] / (2. * l + 1.) \
                            * radius[i] * (radius[i] / radius[j])**(l-1)

                a[i, ilith+1] = 4. * np.pi * g / (2. * l + 1.) \
                    * radius[i] * (radius[i] / r_sigma)**(l-1)

                if tides is True:
                    atides[0, i] = g * mp * radius[i] / rp**3 * (
                        - np.sqrt(5.) / 5. * cp20[0, l, m] +
                        np.sqrt(12./5.) * cp22[0, l, m] / 2.)
                    atides[1, i] = g * mp * radius[i] / rp**3 * (
                        - np.sqrt(5.) / 5. * cp20[1, l, m] +
                        np.sqrt(12./5.) * cp22[1, l, m] / 2.)

            for j in range(1, ilith+1):
                a[ilith+1, j] = 4. * np.pi * g * drho[j] / (2. * l + 1.) \
                    * radius[j] * (radius[j] / r_ref)**(l+1)

            a[ilith+1, ilith+1] = 4. * np.pi * g / (2. * l + 1.)\
                * r_sigma * (r_sigma / r_ref)**(l+1)

            # --- do cosine term ---

            b = np.zeros(ilith+2)
            if l == 2 and m == 0:
                for i in range(1, ilith+1):
                    b[i] = (omega * radius[i])**2 / (3. * np.sqrt(5.))
                    if tides:
                        b[i] += g * mp * radius[i]**2 / rp**3 * \
                            np.sqrt(5.) / 10.
            if l == 2 and m == 2 and tides:
                for i in range(1, ilith+1):
                    b[i] = - g * mp * radius[i]**2 / rp**3 * \
                        np.sqrt(12./5.) / 4.

            # Add contributions from degree 2 relief to degree 4.
            if l == 4 and m <= 2:
                b[1:ilith+1] = b4[0, m, 1:ilith+1]

            for i in range(1, ilith+1):
                b[i] -= gm * cminus[0, l, m] * (radius[i] / r_surface)**l \
                    / r_surface
            b[ilith+1] = gm * potential.coeffs[0, l, m] / r_ref - \
                gm * cplus[0, l, m] * (r_surface / r_ref)**l / r_ref

            # solve the linear equation A h = b
            atemp = a.copy()
            if tides:
                for i in range(1, ilith+1):
                    atemp[i, i] += atides[0, i]

            btemp = b.copy()

            # note that the zero index is not used
            lu, piv, x, info = lapack.dgesv(atemp[1:, 1:], btemp[1:])
            if info != 0:
                raise RuntimeError("lapack.dgesv did not exit properly: {:d}",
                                   info)
            for i in range(1, ilith+1):
                hlm[i].coeffs[0, l, m] = x[i-1]

            # calculate b4 contribution
            if l == 2:
                for i in range(1, ilith+1):
                    if m == 0:
                        b4[0, m, i] = 2. / 3. / np.sqrt(5.) * \
                            omega**2 * radius[i] * \
                            hlm[i].coeffs[0, l, m] * p402020
                    elif m == 1:
                        b4[0, m, i] = 2. / 3. / np.sqrt(5.) * \
                            omega**2 * radius[i] * \
                            hlm[i].coeffs[0, l, m] * p412021
                    elif m == 2:
                        b4[0, m, i] = 2. / 3. / np.sqrt(5.) * \
                            omega**2 * radius[i] * \
                            hlm[i].coeffs[0, l, m] * p422022

            # --- do sine term ---

            b = np.zeros(ilith+2)
            if m != 0:
                # Add contributions from degree 2 relief to degree 4.
                if l == 4 and m <= 2:
                    b[1:ilith+1] = b4[1, m, 1:ilith+1]

                for i in range(1, ilith+1):
                    b[i] -= gm * cminus[1, l, m] * (radius[i] / r_surface)**l \
                            / r_surface
                    b[ilith+1] = gm * potential.coeffs[1, l, m] / r_ref - \
                        gm * cplus[1, l, m] * (r_surface / r_ref)**l / r_ref

                # solve the linear equation A h = b
                atemp = a.copy()
                if tides:
                    for i in range(1, ilith+1):
                        atemp[i, i] += atides[1, i]

                btemp = b.copy()

                # note that the zero index is not used
                lu, piv, x, info = lapack.dgesv(atemp[1:, 1:], btemp[1:])
                if info != 0:
                    raise RuntimeError(
                        "lapack.dgesv did not exit properly: {:d}", info)
                for i in range(1, ilith+1):
                    hlm[i].coeffs[1, l, m] = x[i-1]

                # calculate b4 contribution
                if l == 2:
                    for i in range(1, ilith+1):
                        if m == 1:
                            b4[1, m, i] = 2. / 3. / np.sqrt(5.) * \
                                omega**2 * radius[i] * \
                                hlm[i].coeffs[1, l, m] * p412021
                        elif m == 2:
                            b4[1, m, i] = 2. / 3. / np.sqrt(5.) * \
                                omega**2 * radius[i] * \
                                hlm[i].coeffs[1, l, m] * p422022

    # Calculate potential at r_ref resulting from all interfaces below and
    # including ilith

    coeffs = np.zeros((2, lmax+1, lmax+1))
    for i in range(1, ilith+1):
        for l in range(1, lmax+1):
            coeffs[:, l, :l+1] += hlm[i].coeffs[:, l, :l+1] * 4. * \
                np.pi * drho[i] * radius[i]**2 * (radius[i] / r_ref)**l * \
                g / gm / (2. * l + 1.)

    clm_hydro = pysh.SHGravCoeffs.from_array(coeffs, gm=gm, r0=r_ref,
                                             omega=omega)

    return hlm, clm_hydro, mass_model


def HydrostaticShape(radius, rho, omega, gm, r_ref, rp=None, mp=None,
                     i_clm_hydro=None):
    """
    Calculate the shape of hydrostatic relief in a rotating planet or moon,
    along with the total gravitation potential. For the case of a moon in
    synchronous rotation, optionally include the tidal potential.

    Returns
    -------
    hlm : array of SHCoeffs class instances, size(n+1)
        Array of SHCoeffs class instances of the spherical harmonic
        coefficients of the hydrostatic relief at each interface.
    clm_hydro : SHCoeffs class instance containing the gravitational potential
        resulting from all hydrostatic interfaces. If i_clm_hydro is specified,
        then the potential will include only interfaces beneath index
        i_clm_hydro.
    mass : float
        Total mass of the planet, assuming a spherical shape and the provided
        1D density profile.

    Parameters
    ----------
    radius : ndarray, float, size(n+1)
        Radius of each density interface, where index 0 corresponds to the
        center of the planet and n corresponds to the surface.
    density : ndarray, float, size(n+1)
        Density of layer i between radius[i] and radius[i+1]. The density
        at index 0 is from the center of the planet to radius[1], whereas the
        the density at index n should be zero.
    omega : float
        Angular rotation rate of the planet.
    gm : float
        GM of the planet.
    r_ref : float
        Refernce radius for output potential coefficients.
    rp : float, optional, default = None
        If specified, include the tidal potential acting on a synchronously
        rotating moon, where rp is the average distance between the planet
        and satellite.
    mp : float, optional, default = None
        The mass of the host planet, at a distance rp from the satellite.
    i_clm_hydro : int, optional, default = None
        If specified, calculate the gravitational potential clm_hydro resulting
        from all interfaces below and including the radius index i_clm_hydro.
    """
    tides = False
    if rp is not None:
        if mp is None:
            raise ValueError('When including tides, both rp and mp must be ' +
                             'specified.')
        tides = True
    if mp is not None:
        if rp is None:
            raise ValueError('When including tides, both rp and mp must be ' +
                             'specified.')
        tides = True

    if len(radius) != len(rho):
        raise ValueError('Length of radius and density must be the same.' +
                         'len(radius) = {:d}. len(density) = {:d}.'
                         .format(len(radius), len(rho)))

    n = len(radius) - 1  # index of surface
    lmax = 4

    g = pysh.constants.G.value

    hlm = [pysh.SHCoeffs.from_zeros(lmax) for i in range(n+1)]
    clm_hydro = pysh.SHCoeffs.from_zeros(lmax)

    for i in range(n+1):
        hlm[i].coeffs[0, 0, 0] = radius[i]

    # First determine the spherical harmonic coefficients of (Y20 Ylm)
    # and for tides, (Y22 Ylm). We are only concerned with the coefficient
    # that corresponds to lm, so for each lm, store only the lm component in
    # the array cp20 and cp22.

    sh20 = np.zeros((2, lmax+1, lmax+1))
    sh22 = np.zeros((2, lmax+1, lmax+1))
    sh = np.zeros((2, lmax+1, lmax+1))
    cp20 = np.zeros((2, lmax+1, lmax+1))
    cp22 = np.zeros((2, lmax+1, lmax+1))

    sh20[0, 2, 0] = 1.  # Y20
    sh22[0, 2, 2] = 1.  # Y22

    for l in range(2, lmax+1):
        for m in range(0, l+1):
            sh[0, l, m] = 1.
            coeffs = pysh.expand.SHMultiply(sh20, sh)
            cp20[0, l, m] = coeffs[0, l, m]
            if m != 0:
                cp20[1, l, m] = cp20[0, l, m]
            if l == 2 and m == 0:
                p402020 = coeffs[0, 4, 0]
            if l == 2 and m == 1:
                p412021 = coeffs[0, 4, 1]
            if l == 2 and m == 2:
                p422022 = coeffs[0, 4, 2]

            coeffs = pysh.expand.SHMultiply(sh22, sh)
            cp22[0, l, m] = coeffs[0, l, m]
            sh[0, l, m] = 0.

            if m > 0:
                sh[1, l, m] = 1.
                coeffs = pysh.expand.SHMultiply(sh22, sh)
                cp22[1, l, m] = coeffs[1, l, m]
                sh[1, l, m] = 0.

    # Calculate delta_rho

    drho = np.zeros(n+1)
    mass = np.zeros(n+1)

    for i in range(1, n+1):
        drho[i] = rho[i-1] - rho[i]

    # Calculate matrix A and invert for relief.

    a = np.zeros((n+1, n+1))
    atides = np.zeros((2, n+1))
    b4 = np.zeros((2, 3, n+1))

    # Calculate cumulate mass function
    for i in range(1, n+1):
        if i == 1:
            mass[1] = 4. * np.pi * radius[1]**3 * rho[0] / 3.
        else:
            mass[i] = mass[i-1] + 4. * np.pi * \
                    (radius[i]**3 - radius[i-1]**3) * rho[i-1] / 3.

    mass_model = mass[n]

    for l in range(2, lmax+1, 2):
        for m in range(0, lmax+1):
            for i in range(1, n+1):  # zero index not computed
                for j in range(1, n+1):
                    if i == j:  # cp20 for sin and cosine terms are equal
                        a[i, j] = 4. * np.pi * g * drho[i] * radius[i] / \
                            (2. * l + 1.) - g * mass[i] / radius[i]**2 + \
                            (2./3.) * radius[i] * omega**2 * \
                            (1. - cp20[0, l, m] / np.sqrt(5.0))
                    elif j < i:
                        a[i, j] = 4. * np.pi * g * drho[j] * radius[j] * \
                            (radius[j] / radius[i])**(l+1) / (2. * l + 1.)
                    else:
                        a[i, j] = 4. * np.pi * g * drho[j] * radius[i] * \
                            (radius[i] / radius[j])**(l-1) / (2. * l + 1.)

                if tides is True:
                    atides[0, i] = g * mp * radius[i] / rp**3 * (
                        - np.sqrt(5.) / 5. * cp20[0, l, m] +
                        np.sqrt(12./5.) * cp22[0, l, m] / 2.)
                    atides[1, i] = g * mp * radius[i] / rp**3 * (
                        - np.sqrt(5.) / 5. * cp20[1, l, m] +
                        np.sqrt(12./5.) * cp22[1, l, m] / 2.)

            # --- do cosine term ---

            b = np.zeros(n+1)
            if l == 2 and m == 0:
                for i in range(1, n+1):
                    b[i] = (omega * radius[i])**2 / (3. * np.sqrt(5.))
                    if tides:
                        b[i] += g * mp * radius[i]**2 / rp**3 * \
                            np.sqrt(5.) / 10.
            if l == 2 and m == 2 and tides:
                for i in range(1, n+1):
                    b[i] = - g * mp * radius[i]**2 / rp**3 * \
                        np.sqrt(12./5.) / 4.

            # Add contributions from degree 2 relief to degree 4.
            if l == 4 and m <= 2:
                b[1:n+1] = b4[0, m, 1:n+1]

            # solve the linear equation A h = b
            atemp = a.copy()
            if tides:
                for i in range(1, n+1):
                    atemp[i, i] += atides[0, i]

            btemp = b.copy()

            # note that the zero index is not used
            lu, piv, x, info = lapack.dgesv(atemp[1:, 1:], btemp[1:])
            if info != 0:
                raise RuntimeError(
                    "lapack.dgesv did not exit properly: {:d}", info)
            for i in range(1, n+1):
                hlm[i].coeffs[0, l, m] = x[i-1]

            # calculate b4 contribution
            if l == 2:
                for i in range(1, n+1):
                    if m == 0:
                        b4[0, m, i] = 2. / 3. / np.sqrt(5.) * \
                            omega**2 * radius[i] * \
                            hlm[i].coeffs[0, l, m] * p402020
                    elif m == 1:
                        b4[0, m, i] = 2. / 3. / np.sqrt(5.) * \
                            omega**2 * radius[i] * \
                            hlm[i].coeffs[0, l, m] * p412021
                    elif m == 2:
                        b4[0, m, i] = 2. / 3. / np.sqrt(5.) * \
                            omega**2 * radius[i] * \
                            hlm[i].coeffs[0, l, m] * p422022

            # --- do sine term ---

            b = np.zeros(n+1)
            if m != 0:
                # Add contributions from degree 2 relief to degree 4.
                if l == 4 and m <= 2:
                    b[1:n+1] = b4[1, m, 1:n+1]

                # solve the linear equation A h = b
                atemp = a.copy()
                if tides:
                    for i in range(1, n+1):
                        atemp[i, i] += atides[1, i]

                btemp = b.copy()

                # note that the zero index is not used
                lu, piv, x, info = lapack.dgesv(atemp[1:, 1:], btemp[1:])
                if info != 0:
                    raise RuntimeError(
                        "lapack.dgesv did not exit properly: {:d}", info)
                for i in range(1, n+1):
                    hlm[i].coeffs[1, l, m] = x[i-1]

                # calculate b4 contribution
                if l == 2:
                    for i in range(1, n+1):
                        if m == 1:
                            b4[1, m, i] = 2. / 3. / np.sqrt(5.) * \
                                omega**2 * radius[i] * \
                                hlm[i].coeffs[1, l, m] * p412021
                        elif m == 2:
                            b4[1, m, i] = 2. / 3. / np.sqrt(5.) * \
                                omega**2 * radius[i] * \
                                hlm[i].coeffs[1, l, m] * p422022

    # Calculate potential at r_ref resulting from all interfaces,
    # or only those beneath and including i_clm_hydro.

    coeffs = np.zeros((2, lmax+1, lmax+1))
    if i_clm_hydro is None:
        i_clm_hydro = n

    for i in range(1, i_clm_hydro+1):
        for l in range(2, lmax+1):
            coeffs[:, l, :l+1] += hlm[i].coeffs[:, l, :l+1] * 4. * \
                np.pi * drho[i] * radius[i]**2 * (radius[i] / r_ref)**l * \
                g / gm / (2. * l + 1.)
    coeffs[0, 0, 0] = 1.

    clm_hydro = pysh.SHGravCoeffs.from_array(coeffs, gm=gm, r0=r_ref,
                                             omega=omega)

    return hlm, clm_hydro, mass_model
