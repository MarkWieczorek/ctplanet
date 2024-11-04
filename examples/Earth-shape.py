#!/usr/bin/env python3
"""
Earth-shape

Compute hydrostatic shape of Earth using PREM.
"""
import numpy as np
import pyshtools as pysh

from ctplanet import HydrostaticShape

# ==== MAIN FUNCTION ====


def main():

    model_name = ['PREM.dat']
    spec = 'Data/'
    interior_file = [spec + name for name in model_name]

    r_ref = pysh.constants.Earth.wgs84.a.value
    gm = pysh.constants.Earth.wgs84.gm.value
    mass = gm / pysh.constants.G.value
    omega = pysh.constants.Earth.wgs84.omega.value

    print('Reference radius (km) = {:f}'.format(r_ref / 1.e3))
    print('GM = {:e}'.format(gm))
    print('Mass = {:e}'.format(mass))
    print('Omega = {:e}'.format(omega))

    model = 0

    # --- read 1D reference interior model ---

    print('=== Reading model {:s} ==='.format(model_name[model]))

    with open(interior_file[model], 'r') as f:
        lines = f.readlines()
        print(lines[0].strip())
        data = lines[1].split()
        if float(data[2]) != 1:
            raise RuntimeError('Program not capable of reading polynomial ' +
                               'files')
        num = int(lines[2].split()[0])
        crust_index = int(lines[2].split()[3])
        # mantle_index = int(lines[2].split()[3])
        radius = np.zeros(num)
        rho = np.zeros(num)
        for i in range(0, num):
            data = lines[i+3].split()
            radius[i] = float(data[0])
            rho[i] = float(data[1])

        r0_model = radius[num-1]
        print('Surface radius of model (km) = {:f}'.format(r0_model / 1.e3))

        n = num - 1
        rho[n] = 0.  # the density above the surface is zero
        rho_mantle = rho[crust_index-1]
        print('Mantle density (kg/m3) = {:f}'.format(rho_mantle))

    # --- Compute gravity contribution from hydrostatic density interfaces ---

    hlm_fluid, clm_fluid, mass_model = HydrostaticShape(radius, rho, omega,
                                                        gm, r_ref)

    hlm_surf_unnorm = hlm_fluid[n].convert(normalization='unnorm')
    print('--- Hydrostatic relief of surface ---')
    print('h20 = {:e}\n'.format(hlm_fluid[n].coeffs[0, 2, 0]) +
          'h40 = {:e}'.format(hlm_fluid[n].coeffs[0, 4, 0]))
    print('h20 (unnorm)= {:e}\n'
          .format(hlm_surf_unnorm.coeffs[0, 2, 0]/r0_model) +
          'h40 (unnorm)= {:e}'
          .format(hlm_surf_unnorm.coeffs[0, 4, 0]/r0_model))

    print(mass, mass_model)


# ==== EXECUTE SCRIPT ====


if __name__ == "__main__":
    main()
