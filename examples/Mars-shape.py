#!/usr/bin/env python3
"""
Mars-shape

Calculations and images of the shape of density interfaces of Mars as presented
in Wieczorek et al. (2019).
"""
import os
import numpy as np

import pyshtools as pysh

from ctplanet import HydrostaticShapeLith
from ctplanet import HydrostaticShape
from ctplanet import InertiaTensor_from_shape
from ctplanet import ReadRefModel

# ==== MAIN FUNCTION ====


def main():

    d_lith = 150.e3
    rho_crust = 2900.
    d_sigma = 45.e3
    lmax_hydro = 50

    model_name = ['DWThot', 'DWThotCrust1', 'DWThotCrust1r', 'EH45Tcold',
                  'EH45TcoldCrust1', 'EH45TcoldCrust1r', 'EH45ThotCrust2',
                  'EH45ThotCrust2r', 'LFAK', 'SANAK', 'TAYAK', 'DWAK',
                  'ZG_DW']
    spec = 'Data/Mars-reference-interior-models/Smrekar/'
    interior_file = [spec + name + '.deck' for name in model_name]

    potential = pysh.datasets.Mars.GMM3()
    omega = pysh.constants.Mars.angular_velocity.value
    potential.omega = omega
    mass_mars = np.float64(potential.mass)

    try:
        os.mkdir('figs')
    except Exception:
        pass

    print('Gravity file = {:s}'.format('GMM3'))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}'.format(potential.gm))
    print('Mass = {:e}'.format(potential.mass))
    print('Omega = {:e}'.format(potential.omega))

    lmax_calc = 90
    lmax = 359

    model = 10

    topo = pysh.datasets.Mars.MOLA_shape(lmax=lmax)
    topo.r0 = topo.coeffs[0, 0, 0]
    r_mars = topo.coeffs[0, 0, 0]

    print('Topography file = {:s}'.format('MOLA_shape'))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.r0 / 1.e3))

    # --- Make geoid map, expand in spherical harmonics
    # --- and then remove degree-0 term, for plotting purposes only.
    u0 = potential.gm/potential.r0
    geoid = potential.geoid(u0, r=topo.r0, order=2, lmax_calc=lmax_calc,
                            lmax=lmax)
    geoidsh = geoid.geoid.expand()
    geoidsh.coeffs[0, 0, 0] = 0.
    geoid = geoidsh.expand(grid='DH2')

    # --- read 1D reference interior model ---
    radius, rho, i_crust, i_core, i_lith = ReadRefModel(interior_file[model],
                                                        depth=d_lith)
    # rho_mantle = rho[i_crust-1]
    # rho_core = rho[i_core-1]
    n = len(radius)-1
    # r0_model = radius[n]

    # --- Compute purely hydrostatic relief of all interfaces ---
    print('\n=== Fluid planet ===')
    if False:
        for i in range(1, n+1):
            hlm_fluid, clm_fluid, mass_model = HydrostaticShape(
                radius, rho, omega, potential.gm, potential.r0, i_clm_hydro=i)
            print('i = {:d}, r = {:f}, rho = {:f}, %C20 = {:f}'
                  .format(i, radius[i], rho[i], clm_fluid.coeffs[0, 2, 0] /
                          potential.coeffs[0, 2, 0] * 100))

    hlm_fluid, clm_fluid, mass_model = \
        HydrostaticShape(radius, rho, omega, potential.gm, potential.r0)

    print('--- Hydrostatic relief of surface for a fluid planet ---')
    print('h20 = {:e}\n'.format(hlm_fluid[n].coeffs[0, 2, 0]) +
          'h40 = {:e}'.format(hlm_fluid[n].coeffs[0, 4, 0]))

    hydro_surface = hlm_fluid[n].expand()
    print('Elevation difference between pole and equator (km)'
          ' {:e}'.format(hydro_surface.max()/1.e3 - hydro_surface.min()/1.e3))

    print('--- Hydrostatic relief of core-mantle boundary for a fluid planet '
          '---')
    print('h20 = {:e}\n'.format(hlm_fluid[i_core].coeffs[0, 2, 0]) +
          'h40 = {:e}'.format(hlm_fluid[i_core].coeffs[0, 4, 0]))

    print('Moments of hydrostatic core')
    II, AA, BB, CC, mass, RR, vec = InertiaTensor_from_shape(hlm_fluid, rho,
                                                             i_core,
                                                             quiet=True)
    print('I = ', II)
    print('A, B, C = ', AA, BB, CC)
    print('A, B, C / (mass_mars r0^2) = ', (AA, BB, CC) / mass_mars /
          r_mars**2)
    print('mass of core (kg) = ', mass)
    print('R core (m) = ', RR)

    print('Moments of hydrostatic planet')
    II, AA, BB, CC, mass, RR, vec = InertiaTensor_from_shape(hlm_fluid, rho, n,
                                                             quiet=True)
    print('I = ', II)
    print('A, B, C = ', AA, BB, CC)
    print('A, B, C / (mass_mars r0^2) = ', (AA, BB, CC) / mass_mars /
          r_mars**2)
    print('mass of planet (kg) = ', mass)
    print('R surface (m) = ', RR)

    # --- Compute relief along hydrostatic interfaces with a lithosphere ---
    print('\n=== Planet with a lithosphere ===')
    r_sigma = topo.r0 - d_sigma
    hlm, clm_hydro, mass_model = \
        HydrostaticShapeLith(radius, rho, i_lith, potential, topo, rho_crust,
                             r_sigma, omega, lmax_hydro)
    print('--- Core shape, with lithosphere ---')
    for l in range(0, 5):
        for m in range(0, l+1):
            print(l, m, hlm[i_core].coeffs[0, l, m],
                  hlm[i_core].coeffs[1, l, m], )

    print('Max and min of degree-1 core shape =',
          hlm[i_core].expand(lmax=1).max(),
          hlm[i_core].expand(lmax=1).min())

    II, AA, BB, CC, mass, RR, vec = InertiaTensor_from_shape(hlm, rho, i_core,
                                                             quiet=True)
    print('I = ', II)
    print('A, B, C = ', AA, BB, CC)
    print('A, B, C / (mass_mars r0^2) = ', (AA, BB, CC) / mass_mars /
          r_mars**2)
    print('mass of core (kg) = ', mass)
    print('R core (m) = ', RR)

    # --- Calculate relief, with respect to hydrostatic solution ---
    # --- at i_lith and i_core
    print('\n=== Difference between fluid planet and planet with a '
          'lithosphere ===')
    diff_ilith = hlm[i_lith] - hlm_fluid[i_lith].pad(lmax=lmax_hydro)
    grid_ilith = diff_ilith.expand()
    print('Maximum and minimum difference at i_lith (km) =- ',
          grid_ilith.max()/1.e3, grid_ilith.min()/1.e3)

    diff_icore = hlm[i_core] - hlm_fluid[i_core].pad(lmax=lmax_hydro)
    grid_icore = diff_icore.expand()

    print('Maximum and minimum difference at i_core (km) =- ',
          grid_icore.max()/1.e3, grid_icore.min()/1.e3)

    # ---- Write data to files ---
    print('\n=== Output gridded data for use with GMT ===')
    print('Output grid sizes = ', geoid.nlat, geoid.nlon)

    (geoid/1000).to_file('figs/Mars_geoid.dat')

    diff = (geoid - hlm_fluid[n].pad(lmax=lmax).expand(grid='DH2')
            + radius[n])
    (diff/1000).to_file('figs/Mars_geoid_diff.dat')

    diff = hlm[i_lith].pad(lmax=lmax).expand(grid='DH2') - radius[i_lith]
    (diff/1000).to_file('figs/hydro_ilith.dat')

    diff = hlm[i_core].pad(lmax=lmax).expand(grid='DH2') - radius[i_core]
    (diff/1000).to_file('figs/hydro_icore.dat')

    diff = hlm[i_lith] - hlm_fluid[i_lith].pad(lmax=lmax_hydro)
    (diff.expand(grid='DH2', lmax=lmax)/1000).to_file(
        'figs/hydro_ilith_diff.dat')

    diff = hlm[i_core] - hlm_fluid[i_core].pad(lmax=lmax_hydro)
    (diff.expand(grid='DH2', lmax=lmax)/1000).to_file(
        'figs/hydro_icore_diff.dat')

    # Sensitivity tests
    print('\n=== Sensitivity tests ===')
    r_sigma = topo.r0
    hlm2, clm_hydro2, mass_model2 = \
        HydrostaticShapeLith(radius, rho, i_lith, potential, topo, rho_crust,
                             r_sigma, omega, lmax_hydro)
    print('d_sigma = 0')
    print('Minimum and maximum differences at i_lith (m) = ',
          (hlm2[i_lith] - hlm[i_lith]).expand().min(),
          (hlm2[i_lith] - hlm[i_lith]).expand().max())
    print('Minimum and maximum differences at i_core (m) = ',
          (hlm2[i_core] - hlm[i_core]).expand().min(),
          (hlm2[i_core] - hlm[i_core]).expand().max())

    r_sigma = topo.r0 - 100.e3
    hlm2, clm_hydro2, mass_model2 = \
        HydrostaticShapeLith(radius, rho, i_lith, potential, topo, rho_crust,
                             r_sigma, omega, lmax_hydro)
    print('d_sigma = 100 km')
    print('Minimum and maximum differences at i_lith (m) = ',
          (hlm2[i_lith] - hlm[i_lith]).expand().min(),
          (hlm2[i_lith] - hlm[i_lith]).expand().max())
    print('Minimum and maximum differences at i_core (m) = ',
          (hlm2[i_core] - hlm[i_core]).expand().min(),
          (hlm2[i_core] - hlm[i_core]).expand().max())

    r_sigma = topo.r0 - d_sigma
    hlm2, clm_hydro2, mass_model2 = \
        HydrostaticShapeLith(radius, rho, i_lith, potential, topo, 2500.,
                             r_sigma, omega, lmax_hydro)
    print('rho = 2500 kg / m3')
    print('Minimum and maximum differences at i_lith (m) = ',
          (hlm2[i_lith] - hlm[i_lith]).expand().min(),
          (hlm2[i_lith] - hlm[i_lith]).expand().max())
    print('Minimum and maximum differences at i_core (m) = ',
          (hlm2[i_core] - hlm[i_core]).expand().min(),
          (hlm2[i_core] - hlm[i_core]).expand().max())

    r_sigma = topo.r0 - d_sigma
    hlm2, clm_hydro2, mass_model2 = \
        HydrostaticShapeLith(radius, rho, i_lith, potential, topo, 3300.,
                             r_sigma, omega, lmax_hydro)
    print('rho = 3300 kg / m3')
    print('Minimum and maximum differences at i_lith (m) = ',
          (hlm2[i_lith] - hlm[i_lith]).expand().min(),
          (hlm2[i_lith] - hlm[i_lith]).expand().max())
    print('Minimum and maximum differences at i_core (m) = ',
          (hlm2[i_core] - hlm[i_core]).expand().min(),
          (hlm2[i_core] - hlm[i_core]).expand().max())


# ==== EXECUTE SCRIPT ====


if __name__ == "__main__":
    main()
