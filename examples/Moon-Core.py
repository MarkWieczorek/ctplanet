#!/usr/bin/env python3
"""
Moon-Core

Calculate the hydrostatic shape of the core of the Moon, as presented in
Wieczorek et al. (2019).
"""
import os
import numpy as np
import pyshtools as pysh

from ctplanet import HydrostaticShapeLith
from ctplanet import HydrostaticShape
from ctplanet import InertiaTensor_from_shape

# ==== MAIN FUNCTION ====


def main():

    lmax = 20
    lmax_grid = 719

    omega = pysh.constants.Moon.angular_velocity.value
    rem = pysh.constants.Moon.orbit_semimajor_axis.value
    mass_earth = pysh.constants.Earth.egm2008.mass.value

    cthick = 34.e3  # 43.e3 or 34.0e3
    rho_crust = 2550.

    try:
        os.mkdir('figs')
    except Exception:
        pass

    out_rc_fc = "figs/rc_fc_34_2550.dat"
    out_rc_rhoc = "figs/rc_rhoc_34_2550.dat"
    out_rc_beta = "figs/rc_beta_34_2550.dat"
    sh_core = "figs/core_34.sh"
    sh_core_fluid = "figs/core_fluid_34.sh"
    core_shape_wo_d1 = "figs/core_shape_wo_d1_330_34.dat"
    core_shape = "figs/core_shape_330_34.dat"

    rcore_int = 1.e3
    rcore_start = 250.e3
    rcore_end = 450.e3

    rhocore_start = 5000.
    rhocore_end = 8000.
    rhocore_int = 1.

    potential = pysh.datasets.Moon.GRGM900C()
    topo = pysh.datasets.Moon.LDEM_shape_pa(lmax=900)
    r0 = topo.coeffs[0, 0, 0]

    r_sigma = r0 - cthick

    ismr2 = 0.3927280  # Williams et al. (2014)
    ismr2 = ismr2 * (1738.e3 / r0)**2

    print("Mean planetary radius (km) = {:e}".format(r0 / 1.e3))
    print("Is/MR2 (solid Moon using mean radius) = {:e}".format(ismr2))
    print("Lmax of Gravitational potential = {:d}".format(potential.lmax))
    print("Reference radius of potential model (km) = {:e}"
          .format(potential.r0/1.e3))
    print("GM = {:e}".format(potential.gm))
    mass = potential.gm / pysh.constants.G.value
    print("Mass (kg) = {:e}".format(mass))
    print("Omega = {:e}".format(omega))
    print("Period (days) = {:e}".format(2. * np.pi / omega / 60. / 60. / 24.))
    print("Average crustal thickness (km) = {:e}".format(cthick / 1.e3))
    print("Crustal density (kg/m3) = {:e}".format(rho_crust))

    radius = np.zeros(4)
    radius[0] = 0.
    radius[2] = r0 - cthick
    radius[3] = r0

    rho = np.zeros(4)
    rho[2] = rho_crust

    mass_crust = 4. * np.pi / 3. * rho_crust * (radius[3]**3 - radius[2]**3)

    n = 3
    i_lith = 2
    i_core = 1

    # For each core radius, find rho_mantle and rho_core that fit total mass
    # and moment of inertia

    f_rc_fc = open(out_rc_fc, 'w')
    f_rc_rhoc = open(out_rc_rhoc, 'w')
    f_rc_beta = open(out_rc_beta, 'w')

    for r_core in np.arange(rcore_start, rcore_end + rcore_int, rcore_int,
                            dtype=float):
        radius[1] = r_core
        first = True

        for rho_core in np.arange(rhocore_start, rhocore_end + rhocore_int,
                                  rhocore_int, dtype=float):
            mass_core = 4. * np.pi / 3. * rho_core * r_core**3
            mass_mantle = mass - mass_crust - mass_core
            rho_mantle = mass_mantle * 3. / np.pi / 4. / (radius[2]**3 -
                                                          r_core**3)
            rho[0] = rho_core
            rho[1] = rho_mantle

            if rho_mantle >= rho_core:
                continue

            ismr2_model = moi_solid(radius, rho, n)
            if first is True:
                diff_old = ismr2 - ismr2_model
                first = False
            else:
                diff_new = ismr2 - ismr2_model

                if diff_new * diff_old <= 0.:
                    # interpolate to get the best fitting core density
                    rho_core_final = (rho_core - rhocore_int) - diff_old * \
                        (rhocore_int) / (diff_new - diff_old)
                    rho[0] = rho_core_final
                    mass_core = 4. * np.pi / 3. * rho[0] * r_core**3
                    mass_mantle = mass - mass_crust - mass_core
                    rho_mantle = mass_mantle * 3. / np.pi / 4. / \
                        (radius[2]**3 - r_core**3)
                    rho[1] = rho_mantle

                    hlm, clm_hydro, mass_model = \
                        HydrostaticShapeLith(radius, rho, i_lith,
                                             potential, topo, rho_crust,
                                             r_sigma, omega, lmax,
                                             rp=rem, mp=mass_earth)
                    # set degree-1 terms to zero
                    hlm[1].coeffs[:, 1, :] = 0.

                    a = hlm[1].expand(lat=0., lon=0., lmax_calc=lmax)
                    b = hlm[1].expand(lat=0., lon=90., lmax_calc=lmax)
                    c = hlm[1].expand(lat=90., lon=0., lmax_calc=lmax)

                    f_core = ((a+b)/2. - c) / ((a + b) / 2.)
                    beta_core = (a**2 - b**2) / (a**2 + b**2)

                    print(r_core/1.e3, rho[0], rho[1], f_core, beta_core)
                    f_rc_fc.write('{:e}, {:e}\n'.format(r_core/1.e3, f_core))
                    f_rc_rhoc.write('{:e}, {:e}\n'.format(r_core/1.e3, rho[0]))
                    f_rc_beta.write('{:e}, {:e}\n'
                                    .format(r_core/1.e3, beta_core))

                    if r_core == 330.e3:
                        print("Rcore (km) = {:e}".format(r_core/1.e3))
                        print("A (km) = {:e}".format(a/1.e3))
                        print("B (km) = {:e}".format(b/1.e3))
                        print("C (km) = {:e}".format(c/1.e3))
                        print("rho_core (kg/m3) = {:e}".format(rho[0]))

                        II, AA, BB, CC, mass_model, RR, vec = \
                            InertiaTensor_from_shape(hlm, rho, 1, quiet=False)

                        grid = hlm[i_core].expand(lmax=lmax_grid, grid='DH2')
                        grid.to_file(core_shape_wo_d1)
                        print("Size of output grid = {:d}, {:d}"
                              .format(grid.nlat, grid.nlon))
                        print("Maximum = {:e}\nMinimum = {:e}"
                              .format(grid.max(), grid.min()))

                        hlm, clm_hydro, mass_model = \
                            HydrostaticShapeLith(radius, rho, i_lith,
                                                 potential, topo, rho_crust,
                                                 r_sigma, omega, lmax,
                                                 rp=rem, mp=mass_earth)
                        hlm_fluid, clm_fluid, mass_model = \
                            HydrostaticShape(radius, rho, omega, potential.gm,
                                             potential.r0,
                                             rp=rem, mp=mass_earth)
                        grid = hlm[i_core].expand(lmax=lmax_grid, grid='DH2')
                        grid_fluid = hlm_fluid[i_core].expand(lmax=lmax_grid,
                                                              grid='DH2')
                        a_fluid = hlm_fluid[1].expand(lat=0., lon=0.)
                        b_fluid = hlm_fluid[1].expand(lat=0., lon=90.)
                        c_fluid = hlm_fluid[1].expand(lat=90., lon=0.)

                        grid_fluid_surface = hlm_fluid[3].expand(
                            lmax=lmax_grid, grid='DH2')
                        print("Surface relief: Maximum = {:}, Minimum = {:}"
                              .format(grid_fluid_surface.max(),
                                      grid_fluid_surface.min()))

                        f_core_fluid = (((a_fluid+b_fluid)/2. - c_fluid)
                                        / ((a_fluid + b_fluid) / 2.))
                        beta_core_fluid = ((a_fluid**2 - b_fluid**2) /
                                           (a_fluid**2 + b_fluid**2))
                        print('f_core for a fluid planet = ', f_core_fluid)
                        print('beta_core for a fluid planet = ',
                              beta_core_fluid)

                        diff = grid - grid_fluid
                        print('Maximuim and miniumum core relief for '
                              'a fluid planet (m) = ',
                              grid_fluid.max(), grid_fluid.min())
                        print('Maximuim and miniumum difference with respect '
                              'to a fluid planet (m) = ',
                              diff.max(), diff.min())
                        hlm[i_core].to_file(sh_core)
                        hlm_fluid[i_core].to_file(sh_core_fluid)
                        grid.to_file(core_shape)
                        print("Maximum = {:e}\nMinimum = {:e}"
                              .format(grid.max(), grid.min()))
                        for l in range(0, 4):
                            for m in range(0, l+1):
                                print(l, m, hlm[i_core].coeffs[0, l, m],
                                      hlm[i_core].coeffs[1, l, m])

                diff_old = diff_new


def moi_solid(radius, rho, n):
    """
    Calculate the mean, normalized, moment of inertia of the solid portion
    of the planet.

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
