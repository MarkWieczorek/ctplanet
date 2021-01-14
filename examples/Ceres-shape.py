#!/usr/bin/env python3
"""
Ceres-shape

Calculate the hydrostatic shape of Ceres, as in Wieczorek et al. (2019).
"""
import numpy as np
import pyshtools

from ctplanet import HydrostaticShape

# ==== MAIN FUNCTION ====


def main():

    # Thomas et al. 2005
    mass = 9.395e20
    gm = mass * pyshtools.constants.G.value
    r0 = 476.2e3
    omega = 2 * np.pi / (9.076 * 60 * 60)
    volume = 4 * np.pi / 3 * r0**3

    lmax = 4
    # Model C1
    radius = np.ndarray(2)
    radius[0] = 0.
    radius[1] = r0

    rho = np.ndarray(2)
    rho[0] = mass / volume
    rho[1] = 0

    print(rho[0])
    hlm, clm, mass_model = HydrostaticShape(radius, rho, omega, gm, r0)

    a = hlm[1].expand(lat=0., lon=0., lmax_calc=lmax)
    b = hlm[1].expand(lat=0., lon=90., lmax_calc=lmax)
    c = hlm[1].expand(lat=90., lon=0., lmax_calc=lmax)

    print(a, b, c, a-c)

# ==== EXECUTE SCRIPT ====


if __name__ == "__main__":
    main()
