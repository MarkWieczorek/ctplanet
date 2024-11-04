#!/usr/bin/env python3
"""
Mars-Crust-hydrostatic-tests

Create a crustal thickness map of Mars from gravity and topography and compare
how results change if hydrostatic interfaces are not taken into account.

"""
import os
import numpy as np

import pyshtools as pysh

from ctplanet import pyMoho
from ctplanet import HydrostaticShapeLith
from ctplanet import HydrostaticShape

# ==== MAIN FUNCTION ====


def main():

    model_name = ['DWThot', 'DWThotCrust1', 'DWThotCrust1r', 'EH45Tcold',
                  'EH45TcoldCrust1', 'EH45TcoldCrust1r', 'EH45ThotCrust2',
                  'EH45ThotCrust2r', 'LFAK', 'SANAK', 'TAYAK', 'DWAK',
                  'ZG_DW']
    spec = 'Data/Mars-reference-interior-models/Smrekar/'
    interior_file = [spec + name + '.deck' for name in model_name]

    model = 10

    d_lith = 150.e3

    lmax_calc = 90
    lmax = lmax_calc * 4
    lmax_grid = 719

    potential = pysh.datasets.Mars.GMM3()
    topo = pysh.datasets.Mars.MOLA_shape(lmax=lmax_grid)
    topo.r0 = topo.coeffs[0, 0, 0]

    try:
        os.mkdir('figs')
    except Exception:
        pass

    print('Gravity file = {:s}'.format('GMM3'))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}'.format(potential.gm))
    print('Mass = {:e}'.format(potential.mass))

    print('Topography file = {:s}'.format('MOLA_shape'))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.r0 / 1.e3))

    lat_insight = 4.43
    lon_insight = 135.84

    filter = 1
    half = 50
    nmax = 7
    lmax_hydro = 50
    t0_sigma = 5.  # maximum difference between minimum crustal thickness
    rho_c = 2900.
    t0 = 1.e3  # minimum crustal thickness

    omega = pysh.constants.Mars.angular_velocity.value
    print('Omega = {:e}'.format(omega))

    # --- read 1D reference interior model ---
    with open(interior_file[model], 'r') as f:
        lines = f.readlines()
        print(lines[0].strip())
        data = lines[1].split()
        if float(data[2]) != 1:
            raise RuntimeError('Program not capable of reading polynomial ' +
                               'files')
        num_file = int(lines[2].split()[0])
        crust_index_file = int(lines[2].split()[3])
        core_index_file = int(lines[2].split()[2])
        i_crust_file = crust_index_file - 1
        i_core_file = core_index_file - 1

        radius = np.zeros(num_file)
        rho = np.zeros(num_file)
        num = 0

        for i in range(0, num_file-1):
            data = lines[i+3].split()
            rb = float(data[0])
            rhob = float(data[1])
            data = lines[i+4].split()
            rt = float(data[0])
            rhot = float(data[1])

            if rb == rt:
                if i == i_core_file:
                    i_core = num
                if i == i_crust_file:
                    i_crust = num
            else:
                radius[num] = rb
                rho[num] = (rhot + rhob) / 2.
                num += 1

        radius[num] = rt
        rho[num] = 0.  # the density above the surface is zero
        num += 1
        n = num - 1
        radius = radius[:n+1]
        rho = rho[:n+1]
        r0_model = radius[n]

        print('Surface radius of model (km) = {:f}'.format(r0_model / 1.e3))
        for i in range(0, n+1):
            if radius[i] <= (r0_model - d_lith) and \
                    radius[i+1] > (r0_model - d_lith):
                if radius[i] == (r0_model - d_lith):
                    i_lith = i
                elif (r0_model - d_lith) - radius[i] <= radius[i+1] -\
                        (r0_model - d_lith):
                    i_lith = i
                else:
                    i_lith = i + 1
                break

        rho_mantle = rho[i_crust-1]
        rho_core = rho[i_core-1]
        print('Mantle density (kg/m3) = {:f}'.format(rho_mantle))
        print('Mantle radius (km) = {:f}'.format(radius[i_crust]/1.e3))
        print('Core density (kg/m3) = {:f}'.format(rho_core))
        print('Core radius (km) = {:f}'.format(radius[i_core]/1.e3))

        print('Assumed depth of lithosphere (km) = {:f}'.format(d_lith / 1.e3))
        print('Actual depth of lithosphere in discretized model (km) = {:f}'
              .format((r0_model - radius[i_lith]) / 1.e3))

    # --- Compute gravity contribution from hydrostatic density interfaces ---

    thickave = 44.e3    # initial guess of average crustal thickness
    r_sigma = topo.r0 - thickave

    if True:
        # compute values for a planet that is completely fluid
        hlm_fluid, clm_fluid, mass_model = \
            HydrostaticShape(radius, rho, omega, potential.gm, potential.r0,
                             i_clm_hydro=i_lith)
        print('--- Hydrostatic potential coefficients for a fluid planet ---')
        print('c20 = {:e}\nc40 = {:e}'.format(clm_fluid.coeffs[0, 2, 0],
                                              clm_fluid.coeffs[0, 4, 0]))
        print('--- Hydrostatic relief of surface for a fluid planet ---')
        print('h20 = {:e}\nh40 = {:e}'.format(hlm_fluid[n].coeffs[0, 2, 0],
                                              hlm_fluid[n].coeffs[0, 4, 0]))
        print('% of J2 arising from beneath lithosphere = {:f}'
              .format(clm_fluid.coeffs[0, 2, 0]
                      / potential.coeffs[0, 2, 0] * 100.))

    hlm, clm_hydro, mass_model = \
        HydrostaticShapeLith(radius, rho, i_lith, potential, topo, rho_c,
                             r_sigma, omega, lmax_hydro)

    print('Total mass of model (kg) = {:e}'.format(mass_model))
    print('% of J2 arising from beneath lithosphere = {:f}'
          .format(clm_hydro.coeffs[0, 2, 0]/potential.coeffs[0, 2, 0] * 100.))

    # --- Constant density model without hydrostatic interfaces ---

    print('-- Constant density model --\nrho_c = {:f}'.format(rho_c))

    tmin = 1.e9
    thickave = 44.e3    # initial guess of average crustal thickness

    while abs(tmin - t0) > t0_sigma:
        # iterate to fit assumed minimum crustal thickness

        moho = pyMoho(potential, topo, lmax, rho_c, rho_mantle,
                      thickave, filter_type=filter, half=half,
                      lmax_calc=lmax_calc, nmax=nmax, quiet=True)

        thick_grid = (topo.pad(lmax_grid)
                      - moho.pad(lmax_grid)).expand(grid='DH2')
        print('Average crustal thickness (km) = {:f}'.format(thickave / 1.e3))
        print('Crustal thickness at InSight landing sites (km) = {:f}'
              .format((topo.pad(lmax) - moho.pad(lmax))
                      .expand(lat=lat_insight, lon=lon_insight) / 1.e3))
        tmin = thick_grid.min()
        tmax = thick_grid.max()
        print('Minimum thickness (km) = {:e}'.format(tmin / 1.e3))
        print('Maximum thickness (km) = {:e}'.format(tmax / 1.e3))
        thickave += t0 - tmin

    (thick_grid/1.e3).plot(show=False, colorbar='bottom',
                           cb_label='Crustal thickness, km',
                           fname='figs/Thick-Mars-without-hydro.png')

    print('Thickness at north pole (km) = ', thick_grid.data[0, 0] / 1.e3)

    # --- Constant density model with hydrostatic interfaces beneath
    # --- lithosphere

    potential.coeffs[:, :lmax_hydro+1, :lmax_hydro+1] -= \
        clm_hydro.coeffs[:, :lmax_hydro+1, :lmax_hydro+1]

    print('-- Constant density model with hydrostatic interfaces beneath '
          'lithosphere--\nrho_c = {:f}'.format(rho_c))

    tmin = 1.e9
    thickave = 44.e3    # initial guess of average crustal thickness

    while abs(tmin - t0) > t0_sigma:
        # iterate to fit assumed minimum crustal thickness

        moho = pyMoho(potential, topo, lmax, rho_c, rho_mantle,
                      thickave, filter_type=filter, half=half,
                      lmax_calc=lmax_calc, nmax=nmax, quiet=True)

        thick2_grid = (topo.pad(lmax_grid)
                       - moho.pad(lmax_grid)).expand(grid='DH2')
        print('Average crustal thickness (km) = {:f}'.format(thickave / 1.e3))
        print('Crustal thickness at InSight landing sites (km) = {:f}'
              .format((topo.pad(lmax) - moho.pad(lmax))
                      .expand(lat=lat_insight, lon=lon_insight) / 1.e3))
        tmin = thick2_grid.min()
        tmax = thick2_grid.max()
        print('Minimum thickness (km) = {:e}'.format(tmin / 1.e3))
        print('Maximum thickness (km) = {:e}'.format(tmax / 1.e3))
        thickave += t0 - tmin

    (thick2_grid/1.e3).plot(show=False, colorbar='bottom',
                            cb_label='Crustal thickness, km',
                            fname='figs/Thick-Mars-with-hydro-lith.png')
    (thick2_grid/1.e3 - thick_grid/1.e3).plot(
        show=False, colorbar='bottom', cb_label='Crustal thickness, km',
        fname='figs/Thick-Mars-diff-hydro-lith.png')
    print('Thickness at north pole (km) = ', thick2_grid.data[0, 0] / 1.e3)
    min = (thick2_grid/1.e3 - thick_grid/1.e3).min()
    max = (thick2_grid/1.e3 - thick_grid/1.e3).max()
    print('Minimum and maximum difference (km) = ', min, max)
    (thick2_grid/1.e3).to_file('figs/Thick-Mars-hydro-lith.dat')
    (thick2_grid/1.e3 - thick_grid/1.e3).to_file(
        'figs/Thick-Mars-diff-hydro-lith.dat')

    # --- Constant density model with hydrostatic interfaces of fluid planet

    potential.coeffs[:, :lmax_hydro+1, :lmax_hydro+1] += \
        clm_hydro.coeffs[:, :lmax_hydro+1, :lmax_hydro+1]
    potential.coeffs[:, :5, :5] -= clm_fluid.coeffs[:, :5, :5]

    print('-- Constant density model with hydrostatic interfaces of a fluid '
          'planet--\nrho_c = {:f}'.format(rho_c))

    tmin = 1.e9
    thickave = 44.e3    # initial guess of average crustal thickness

    while abs(tmin - t0) > t0_sigma:
        # iterate to fit assumed minimum crustal thickness

        moho = pyMoho(potential, topo, lmax, rho_c, rho_mantle,
                      thickave, filter_type=filter, half=half,
                      lmax_calc=lmax_calc, nmax=nmax, quiet=True)

        thick3_grid = (topo.pad(lmax_grid)
                       - moho.pad(lmax_grid)).expand(grid='DH2')
        print('Average crustal thickness (km) = {:f}'.format(thickave / 1.e3))
        print('Crustal thickness at InSight landing sites (km) = {:f}'
              .format((topo.pad(lmax) - moho.pad(lmax))
                      .expand(lat=lat_insight, lon=lon_insight) / 1.e3))
        tmin = thick3_grid.min()
        tmax = thick3_grid.max()
        print('Minimum thickness (km) = {:e}'.format(tmin / 1.e3))
        print('Maximum thickness (km) = {:e}'.format(tmax / 1.e3))
        thickave += t0 - tmin

    (thick3_grid/1.e3).plot(show=False, colorbar='bottom',
                            cb_label='Crustal thickness, km',
                            fname='figs/Thick-Mars-with-fluid.png')
    print('Thickness at north pole (km) = ', thick3_grid.data[0, 0] / 1.e3)
    (thick2_grid/1.e3 - thick3_grid/1.e3).plot(
        show=False, colorbar='bottom', cb_label='Crustal thickness, km',
        fname='figs/Thick-Mars-diff-lith-fluid.png')
    min = (thick2_grid/1.e3 - thick3_grid/1.e3).min()
    max = (thick2_grid/1.e3 - thick3_grid/1.e3).max()
    print('Minimum and maximum difference (km) = ', min, max)
    (thick2_grid/1.e3 - thick3_grid/1.e3).to_file(
        'figs/Thick-Mars-diff-lith-fluid.dat')

# ==== EXECUTE SCRIPT ====


if __name__ == "__main__":
    main()
