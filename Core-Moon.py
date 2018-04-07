#!/usr/bin/env python3
"""
    Calculate the hydrostatic shape of the core of the Moon.

"""
import numpy as np

import pyshtools

from HydrostaticShapeLith import HydrostaticFlatteningLith

# ==== MAIN FUNCTION ====


def main():

    lmax_grid = 359
    lmax = 6

    omega = float(pyshtools.constant.omega_moon)
    rem = float(pyshtools.constant.a_moon)
    mass_earth = float(pyshtools.constant.mass_earth)
    r0 = float(pyshtools.constant.r_moon)

    cthick = 34.e3  # 43.e3 or 34.0e3
    rho_crust = 2550.

    out_rc_fc = "rc_fc_34_2550.dat"
    out_rc_rhoc = "rc_rhoc_34_2550.dat"
    out_rc_beta = "rc_beta_34_2550.dat"
    sh_core = "core_34.sh"
    core_shape = "core_shape_330_34.dat"

    rcore_int = 1.e3
    rcore_start = 250.e3
    rcore_end = 450.e3

    rhocore_start = 5000.
    rhocore_end = 8000.
    rhocore_int = 1.

    ismr2 = 0.3927280  # Williams et al. (2014)
    ismr2 = ismr2 * (1738.e3 / r0)**2

    pot_file = "/Users/lunokhod/Moon/GRAIL/GravityModels/" + \
        "JGGRAIL_900C11A_SHA.TAB"

    coeffs, lmaxp, header = pyshtools.shio.shread(pot_file, lmax=10,
                                                  header=True)
    r0_pot = float(header[0])*1.e3
    gm = float(header[1])*1.e9

    potential = pyshtools.SHCoeffs.from_array(coeffs)
    potential.r_ref = r0_pot
    potential.gm = gm

    print("Mean planetary radius (km) = {:e}".format(r0 / 1.e3))
    print("Is/MR2 (solid Moon using mean radius) = {:e}".format(ismr2))
    print("Lmax of Gravitational potential = {:d}".format(lmaxp))
    print("Reference radius of potential model (km) = {:e}"
          .format(r0_pot/1.e3))
    print("GM = {:e}".format(gm))
    mass = gm / float(pyshtools.constant.grav_constant)
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
                        HydrostaticFlatteningLith(radius, rho, i_lith,
                                                  potential, omega, lmax,
                                                  finiteamplitude=False,
                                                  rp=rem, mp=mass_earth)

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
                        hlm[i_core].to_file(sh_core)

                        print("Rcore (km) = {:e}".format(r_core/1.e3))
                        print("A (km) = {:e}".format(a/1.e3))
                        print("B (km) = {:e}".format(b/1.e3))
                        print("C (km) = {:e}".format(c/1.e3))
                        print("rho_core (kg/m3) = {:e}".format(rho[0]))

                        grid = hlm[i_core].expand(lmax=lmax_grid, grid='DH2')
                        grid.to_file(core_shape)
                        print("Size of output grid = {:d}, {:d}"
                              .format(grid.nlat, grid.nlon))
                        print("Maximum = {:e}\nMinimum = {:e}"
                              .format(grid.data.max(), grid.data.min()))

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
