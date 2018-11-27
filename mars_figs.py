#!/usr/bin/env python3
"""
Create images related to Mars for the Hydrostatic flattening paper with a
lithosphere
"""
import numpy as np

import pyshtools

import pyMoho
from Hydrostatic import HydrostaticShapeLith
from Hydrostatic import HydrostaticShape
from InertiaTensor import InertiaTensor

# ==== MAIN FUNCTION ====


def main():

    d_lith = 150.e3
    rho_crust = 2900.
    d_sigma = 45.e3
    lmax_hydro = 90

    gravfile = 'Data/gmm3_120_sha.tab'
    topofile = 'Data/MarsTopo719.shape'
    densityfile = 'Data/dichotomy_359.sh'

    model_name = ['DWThot', 'DWThotCrust1', 'DWThotCrust1r', 'EH45Tcold',
                  'EH45TcoldCrust1', 'EH45TcoldCrust1r', 'EH45ThotCrust2',
                  'EH45ThotCrust2r', 'LFAK', 'SANAK', 'TAYAK', 'DWAK',
                  'ZG_DW']
    spec = 'Data/Mars-reference-interior-models/Smrekar/'
    interior_file = [spec + name + '.deck' for name in model_name]

    potential = pyshtools.SHGravCoeffs.from_file(gravfile, header_units='km')
    omega = pyshtools.constant.omega_mars.value
    potential.omega = omega

    print('Gravity file = {:s}'.format(gravfile))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}'.format(potential.gm))
    print('Mass = {:e}'.format(potential.mass))
    print('Omega = {:e}'.format(potential.omega))

    lmax_calc = 90
    lmax = 359

    model = 10

    topo = pyshtools.SHCoeffs.from_file(topofile, lmax=lmax)
    topo.r0 = topo.coeffs[0, 0, 0]

    print('Topography file = {:s}'.format(topofile))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.r0 / 1.e3))

    # --- Make geoid map, expand in spherical harmonics
    # --- and then remove degree-0 term
    u0 = potential.gm/potential.r0
    geoid = potential.geoid(u0,
                            r=topo.r0, order=2,
                            lmax_calc=lmax_calc, lmax=lmax)
    geoidsh = geoid.geoid.expand()
    geoidsh.coeffs[0, 0, 0] = 0.
    geoid = geoidsh.expand(grid='DH2')

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

    # --- Compute purely hydrostatic relief of all interfaces ---
    if False:
        for i in range(1, n+1):
            hlm_fluid, clm_fluid, mass_model = HydrostaticShape(
                radius, rho, omega, potential.gm, potential.r0, i_clm_hydro=i)
            print('i = {:d}, r = {:f}, rho = {:f}, %C20 = {:f}'
                  .format(i, radius[i], rho[i], clm_fluid.coeffs[0, 2, 0] /
                          potential.coeffs[0, 2, 0] * 100))

    hlm_fluid, clm_fluid, mass_model = \
        HydrostaticShape(radius, rho, omega, potential.gm, potential.r0)

    print('--- Hydrostatic relief of surface ---')
    print('h20 = {:e}\n'.format(hlm_fluid[n].coeffs[0, 2, 0]) +
          'h40 = {:e}'.format(hlm_fluid[n].coeffs[0, 4, 0]))

    hydro_surface = hlm_fluid[n].expand()
    print('Elevation difference between pole and equator (km)'
          ' {:e}'.format(hydro_surface.max()/1.e3 - hydro_surface.min()/1.e3))

    II, AA, BB, CC, mass, RR, vec = InertiaTensor(hlm_fluid, rho, i_core,
                                                  quiet=False)
    print('Moments of hydrostatic core')
    print('I = ', II)
    print('A, B, C = ', AA, BB, CC)
    print('mass of core (kg) = ', mass)
    print('R core (m) = ', RR)

    II, AA, BB, CC, mass, RR, vec = InertiaTensor(hlm_fluid, rho, n,
                                                  quiet=False)
    print('Moments of hydrostatic planet')
    print('I = ', II)
    print('A, B, C = ', AA, BB, CC)
    print('mass of planet (kg) = ', mass)
    print('R surface (m) = ', RR)

    # --- Compute relief along hydrostatic interfaces with a lithosphere ---
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

    # Calculate moments of inertial of the core, assuming that A and B are in
    # the plane of the equator, and that C is aligned with the rotation axis.
    II, AA, BB, CC, mass, RR, vec = InertiaTensor(hlm, rho, i_core,
                                                  quiet=False)
    print('I = ', II)
    print('A, B, C = ', AA, BB, CC)
    print('mass of core (kg) = ', mass)
    print('R core (m) = ', RR)

    # --- Calculate relief, with respect to hydrostatic solution ---
    # --- at i_lith and i_core
    diff_ilith = hlm[i_lith] - hlm_fluid[i_lith].pad(lmax=lmax_hydro)
    grid_ilith = diff_ilith.expand()
    print('Maximum and minimum difference at i_lith (km) =- ',
          grid_ilith.max()/1.e3, grid_ilith.min()/1.e3)

    diff_icore = hlm[i_core] - hlm_fluid[i_core].pad(lmax=lmax_hydro)
    grid_icore = diff_icore.expand()

    print('Maximum and minimum difference at i_core (km) =- ',
          grid_icore.max()/1.e3, grid_icore.min()/1.e3)

    # ---- Write data to files ---
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


def moi(radius, rho, n):
    """
    Calculate the mean, normalized, moment of inertia up to index n.

    The radius and density are discretized into shells as in the hydrostatic
    flattening routines:
        radius[0] = 0
        radius[1] = radius of core
        radius[n] = surface
        rho[i] = density from radius[i] to radius[i+1]
    """
    moi_solid = 0.
    mass = 4. * np.pi / 3. * rho[0] * radius[1]**3

    for i in range(2, n+1):
        mass += 4. * np.pi / 3. * rho[i-1] * (radius[i]**3 - radius[i-1]**3)

        moi_solid += 8. * np.pi / 15. * rho[i-1] * (radius[i]**5 -
                                                    radius[i-1]**5)

    return moi_solid / mass / radius[n]**2

# ==== EXECUTE SCRIPT ====


if __name__ == "__main__":
    main()
