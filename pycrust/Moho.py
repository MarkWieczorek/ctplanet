'''
Functions for calculating the relief along the crust-mantle interface of
a planet using gravity and topography.
'''
import numpy as np

import pyshtools as pysh


# ==== pyMoho ====


def pyMoho(pot, topo, lmax, rho_c, rho_m, thickave, filter_type=0, half=None,
           nmax=8, delta_max=5., lmax_calc=None, correction=None, quiet=False):
    """
    Calculate the relief along the crust-mantle interface assuming a
    constant density crust and mantle.

    Returns
    -------
    moho : SHCoeffs class instance
        The radius of the crust-mantle interface.

    Parameters
    ----------
    pot : SHGravCoeffs class instance
        Gravitational potential spherical harmonic coefficients.
    topo : SHCoeffs class instance
        Spherical harmonic coefficients of the surface relief.
    lmax : int
        Maximum spherical harmonic degree of the function, which determines the
        sampling interval of the internally computed grids.
    rho_c : float
        Crustal density in kg / m^3.
    rho_m : float
        Mantle density in kg / m^3.
    thickave : float
        Average thickness of the crust in meters.
    filter_type : int, optional, default = 0
        0 = no filtering, 1 = minimum amplitude filter, 2 = minimum
        curvature filter.
    half : float, optional, default = None
        The spherical harmonic degree where the filter is equal to 0.5. This
        must be set when filter_type is 1 or 2.
    nmax : int, optional, default = 8
        The maximum order used in the Taylor-series expansion when calculating
        the potential coefficients.
    delta_max : float, optional, default = 5.0
        The algorithm will continue to iterate until the maximum difference in
        relief between solutions is less than this value (in meters).
    lmax_calc : int
        Maximum spherical harmonic degree when evalulating the functions.
    correction : SHGravCoeffs class instance, optional, default = None
        If present, these coefficients will be added to the Bouguer correction
        (subtracted from the Bouguer anomaly) before performing the inversion.
        This could be used to account for the pre-computed gravitational
        attraction of the polar caps of Mars, which have a different density
        than the crust.
    quiet : boolean, optional, default = False
        If True, suppress printing output during the iterations.
    """
    if (filter_type == 1 or filter_type == 2) and half is None:
        raise ValueError("half must be set when filter_type is either 1 or 2.")

    if lmax_calc is None:
        lmax_calc = lmax
    d = topo.coeffs[0, 0, 0] - thickave

    mass = pot.mass

    topo_grid = topo.expand(grid='DH2', lmax=lmax, extend=False)

    if quiet is False:
        print("Maximum radius (km) = {:f}".format(topo_grid.data.max() / 1.e3))
        print("Minimum radius (km) = {:f}".format(topo_grid.data.min() / 1.e3))

    bc, r0 = pysh.gravmag.CilmPlusDH(topo_grid.data, nmax, mass,
                                     rho_c, lmax=lmax_calc)
    if correction is not None:
        bc += correction.change_ref(r0=r0).to_array(lmax=lmax_calc)

    pot2 = pot.change_ref(r0=r0)
    ba = pot2.to_array(lmax=lmax_calc, errors=False) - bc

    moho = pysh.SHCoeffs.from_zeros(lmax=lmax_calc)
    moho.coeffs[0, 0, 0] = d

    for l in range(1, lmax_calc + 1):
        if filter_type == 0:
            moho.coeffs[:, l, :l + 1] = ba[:, l, :l + 1] * mass * \
                (2 * l + 1) * ((r0 / d)**l) / \
                (4. * np.pi * (rho_m - rho_c) * d**2)
        elif filter_type == 1:
            moho.coeffs[:, l, :l + 1] = pysh.gravmag.DownContFilterMA(
                l, half, r0, d) * ba[:, l, :l + 1] * mass * \
                (2 * l + 1) * ((r0 / d)**l) / \
                (4. * np.pi * (rho_m - rho_c) * d**2)
        else:
            moho.coeffs[:, l, :l + 1] = pysh.gravmag.DownContFilterMC(
                l, half, r0, d) * ba[:, l, :l + 1] * mass * \
                (2 * l + 1) * ((r0 / d)**l) / \
                (4. * np.pi * (rho_m - rho_c) * d**2)

    moho_grid3 = moho.expand(grid='DH2', lmax=lmax, lmax_calc=lmax_calc,
                             extend=False)

    temp_grid = topo_grid - moho_grid3
    if quiet is False:
        print('Maximum Crustal thickness (km) = {:e}'.format(
            temp_grid.data.max() / 1.e3))
        print('Minimum Crustal thickness (km) = {:e}'.format(
            temp_grid.data.min() / 1.e3))

    moho.coeffs = pysh.gravmag.BAtoHilmDH(ba, moho_grid3.data, nmax,
                                          mass, r0, (rho_m - rho_c),
                                          lmax=lmax,
                                          filter_type=filter_type,
                                          filter_deg=half,
                                          lmax_calc=lmax_calc)

    moho_grid2 = moho.expand(grid='DH2', lmax=lmax, lmax_calc=lmax_calc,
                             extend=False)
    temp_grid = topo_grid - moho_grid2

    if quiet is False:
        print('Delta (km) = {:e}'.format(abs(moho_grid3.data -
                                             moho_grid2.data).max() / 1.e3))
        print('Maximum Crustal thickness (km) = {:f}'
              .format(temp_grid.data.max() / 1.e3))
        print('Minimum Crustal thickness (km) = {:f}'
              .format(temp_grid.data.min() / 1.e3))

    iter = 0
    delta = 1.0e9

    while delta > delta_max:
        iter += 1

        if quiet is False:
            print('Iteration {:d}'.format(iter))

        moho_grid = (moho_grid2 + moho_grid3) / 2.
        temp_grid = topo_grid - moho_grid

        if quiet is False:
            print("Delta (km) = {:e}".format(
                abs(moho_grid.data - moho_grid2.data).max() / 1.e3))
            print('Maximum Crustal thickness (km) = {:f}'.format(
                temp_grid.data.max() / 1.e3))
            print('Minimum Crustal thickness (km) = {:f}'.format(
                temp_grid.data.min() / 1.e3))

        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        iter += 1

        if quiet is False:
            print('Iteration {:d}'.format(iter))

        moho.coeffs = pysh.gravmag.BAtoHilmDH(ba, moho_grid2.data, nmax,
                                              mass, r0,
                                              (rho_m - rho_c), lmax=lmax,
                                              filter_type=filter_type,
                                              filter_deg=half,
                                              lmax_calc=lmax_calc)

        moho_grid = moho.expand(grid='DH2', lmax=lmax, lmax_calc=lmax_calc,
                                extend=False)

        delta = abs(moho_grid.data - moho_grid2.data).max()
        temp_grid = topo_grid - moho_grid

        if quiet is False:
            print('Delta (km) = {:e}'.format(delta / 1.e3))
            print('Maximum Crustal thickness (km) = {:f}'.format(
                temp_grid.data.max() / 1.e3))
            print('Minimum Crustal thickness (km) = {:f}'.format(
                temp_grid.data.min() / 1.e3))

        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        if abs(temp_grid.data).max() > 500.e3:
            print('Not converging')
            exit(1)

    return moho


# ==== pyMohoRho ====


def pyMohoRho(pot, topo, density, porosity, lmax, rho_m, thickave,
              filter_type=0, half=None, nmax=8, delta_max=5., lmax_calc=None,
              correction=None, quiet=False):
    """
    Calculate the relief along the crust-mantle interface assuming a
    constant density mantle and a laterally varying crustal density.

    Returns
    -------
    moho : SHCoeffs class instance
        The radius of the crust-mantle interface.

    Parameters
    ----------
    pot : SHGravCoeffs class instance
        Gravitational potential spherical harmonic coefficients.
    topo : SHCoeffs class instance
        Spherical harmonic coefficients of the surface relief.
    density : SHCoeffs class instance
        Spherical harmonic coefficients of the crustal grain density.
    porosity : float
        Crustal porosity (from 0 to 1).
    lmax : int
        Maximum spherical harmonic degree of the function, which determines the
        sampling interval of the internally computed grids.
    rho_m : float
        Mantle density in kg / m^3.
    thickave : float
        Average thickness of the crust in meters.
    filter_type : int, optional, default = 0
        0 = no filtering, 1 = minimum amplitude filter, 2 = minimum
        curvature filter.
    half : float, optional, default = None
        The spherical harmonic degree where the filter is equal to 0.5. This
        must be set when filter_type is 1 or 2.
    nmax : int, optional, default = 8
        The maximum order used in the Taylor-series expansion when calculating
        the potential coefficients.
    delta_max : float, optional, default = 5.0
        The algorithm will continue to iterate until the maximum difference in
        relief between solutions is less than this value (in meters).
    lmax_calc : int
        Maximum spherical harmonic degree when evalulating the functions.
    correction : SHGravCoeffs class instance, optional, default = None
        If present, these coefficients will be added to the Bouguer correction
        (subtracted from the Bouguer anomaly) before performing the inversion.
        This could be used to account for the pre-computed gravitational
        attraction of the polar caps of Mars, which have a different density
        than the crust.
    quiet : boolean, optional, default = False
        If True, suppress printing output during the iterations.
    """
    if (filter_type == 1 or filter_type == 2) and half is None:
        raise ValueError("half must be set when filter_type is either 1 or 2.")

    if lmax_calc is None:
        lmax_calc = lmax

    d = topo.coeffs[0, 0, 0] - thickave
    rho_crust_ave = density.coeffs[0, 0, 0] * (1. - porosity)

    mass = pot.mass

    topo_grid = topo.expand(grid='DH2', lmax=lmax, extend=False)
    density_grid = density.expand(grid='DH2', lmax=lmax, extend=False)

    if quiet is False:
        print("Maximum radius (km) = {:f}".format(topo_grid.data.max() / 1.e3))
        print("Minimum radius (km) = {:f}".format(topo_grid.data.min() / 1.e3))
        print("Maximum density (kg/m3) = {:f}".format(
            density_grid.data.max() / 1.e3))
        print("Minimum desntiy (kg/m3) = {:f}".format(
            density_grid.data.min() / 1.e3))

    bc, r0 = pysh.gravmag.CilmPlusRhoHDH(
        topo_grid.data, nmax, mass, density_grid.data * (1. - porosity),
        lmax=lmax_calc)
    if correction is not None:
        bc += correction.change_ref(r0=r0).to_array(lmax=lmax_calc)

    pot2 = pot.change_ref(r0=r0)
    ba = pot2.to_array(lmax=lmax_calc, errors=False) - bc

    # next subtract lateral variations in the crust without reflief
    for l in range(1, lmax_calc + 1):
        ba[:, l, :l + 1] = ba[:, l, :l + 1] \
                           - 4. * np.pi * density.coeffs[:, l, :l + 1] \
                           * (1. - porosity) \
                           * (r0**3 - (d**3)*(d/r0)**l) \
                           / (2 * l + 1) / (l + 3) / mass

    moho = pysh.SHCoeffs.from_zeros(lmax=lmax_calc)
    moho.coeffs[0, 0, 0] = d

    for l in range(1, lmax_calc + 1):
        if filter_type == 0:
            moho.coeffs[:, l, :l + 1] = ba[:, l, :l + 1] * mass * \
                (2 * l + 1) * ((r0 / d)**l) / \
                (4. * np.pi * (rho_m - rho_crust_ave) * d**2)
        elif filter_type == 1:
            moho.coeffs[:, l, :l + 1] = pysh.gravmag.DownContFilterMA(
                l, half, r0, d) * ba[:, l, :l + 1] * mass * \
                (2 * l + 1) * ((r0 / d)**l) / \
                (4. * np.pi * (rho_m - rho_crust_ave) * d**2)
        else:
            moho.coeffs[:, l, :l + 1] = pysh.gravmag.DownContFilterMC(
                l, half, r0, d) * ba[:, l, :l + 1] * mass * \
                (2 * l + 1) * ((r0 / d)**l) / \
                (4.0 * np.pi * (rho_m - rho_crust_ave) * d**2)

    moho_grid3 = moho.expand(grid='DH2', lmax=lmax, lmax_calc=lmax_calc,
                             extend=False)
    drho_grid = rho_m - density_grid * (1. - porosity)
    temp_grid = topo_grid - moho_grid3

    if quiet is False:
        print('Maximum Crustal thickness (km) = {:f}'.format(
            temp_grid.data.max() / 1.e3))
        print('Minimum Crustal thickness (km) = {:f}'.format(
            temp_grid.data.min() / 1.e3))

    moho.coeffs = pysh.gravmag.BAtoHilmRhoHDH(
        ba, moho_grid3.data, drho_grid.data, nmax, mass, r0,
        lmax=lmax, filter_type=filter_type, filter_deg=half,
        lmax_calc=lmax_calc)

    moho_grid2 = moho.expand(grid='DH2', lmax=lmax, lmax_calc=lmax_calc,
                             extend=False)
    temp_grid = topo_grid - moho_grid2

    if quiet is False:
        print('Delta (km) = {:e}'.format(abs(moho_grid3.data -
                                             moho_grid2.data).max() / 1.e3))
        print('Maximum Crustal thickness (km) = {:f}'
              .format(temp_grid.data.max() / 1.e3))
        print('Minimum Crustal thickness (km) = {:f}'
              .format(temp_grid.data.min() / 1.e3))

    iter = 0
    delta = 1.0e9

    while delta > delta_max:
        iter += 1

        if quiet is False:
            print('Iteration {:d}'.format(iter))

        moho_grid = (moho_grid2 + moho_grid3) / 2.
        temp_grid = topo_grid - moho_grid

        if quiet is False:
            print("Delta (km) = {:e}".format(
                abs(moho_grid.data - moho_grid2.data).max() / 1.e3))
            print('Maximum Crustal thickness (km) = {:e}'.format(
                temp_grid.data.max() / 1.e3))
            print('Minimum Crustal thickness (km) = {:e}'.format(
                temp_grid.data.min() / 1.e3))

        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        iter += 1

        if quiet is False:
            print('Iteration {:d}'.format(iter))

        moho.coeffs = pysh.gravmag.BAtoHilmRhoHDH(
            ba, moho_grid2.data, drho_grid.data, nmax, mass, r0,
            lmax=lmax, filter_type=filter_type, filter_deg=half,
            lmax_calc=lmax_calc)

        moho_grid = moho.expand(grid='DH2', lmax=lmax, lmax_calc=lmax_calc,
                                extend=False)

        delta = abs(moho_grid.data - moho_grid2.data).max()
        temp_grid = topo_grid - moho_grid

        if quiet is False:
            print('Delta (km) = {:e}'.format(delta / 1.e3))
            print('Maximum Crustal thickness (km) = {:f}'.format(
                temp_grid.data.max() / 1.e3))
            print('Minimum Crustal thickness (km) = {:f}'.format(
                temp_grid.data.min() / 1.e3))

        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        if abs(temp_grid.data).max() > 500.e3:
            print('Not converging')
            exit(1)

    return moho
