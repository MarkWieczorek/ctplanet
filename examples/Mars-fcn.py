#!/usr/bin/env python3
"""
Mars-fcn

Compute the free core nutation period of Mars.
"""
import numpy as np
import pyshtools as pysh

from ctplanet import InertiaTensor_from_C


# ==== MAIN FUNCTION ====


def main():

    lmax = 2

    potential = pysh.datasets.Mars.GMM3(lmax=lmax)

    mass_mars = potential.mass
    r0_pot = potential.r0
    omega = pysh.constants.Mars.angular_velocity.value

    print('Gravity file = {:s}'.format('GMM3'))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}'.format(potential.gm))
    print('G, m3 / (kg s2) = {:e}'.format(pysh.constants.G.value))
    print('Mass = {:e}'.format(potential.mass))
    print('Omega = {:e}'.format(omega))

    topo = pysh.datasets.Mars.MOLA_shape(lmax=lmax)
    r_mars = topo.coeffs[0, 0, 0]

    print('Topography file = {:s}'.format('MOLA_shape'))
    print('Mean planetary radius (km) = {:f}\n'.format(r_mars / 1.e3))

    # Compute polar moment using precession rate and other constants
    # defined in Konopliv et al. 2016 (and 2011).

    J2 = - potential.to_array(normalization='unnorm', errors=False)[0, 2, 0]
    print('J2 = {:e}'.format(J2))
    phidot = -7608.3  # +- 2.1 mas / yr (Konopliv et al. 2016)
    phidotp = -0.2  # planetary torque correction, from Konopliv et al. 2011
    phidotg = 6.7  # relativistic correction, from Konopliv et al. 2011
    phi0dot = phidot - phidotp - phidotg
    print('Phi0-dot, mas/yr = {:f}'.format(phi0dot))
    phi0dot *= np.pi / 180 / 1000 / 60 / 60 / 365.25 / 24 / 60 / 60

    e_mars = 0.09341  # Konopliv et al. 2011
    print('e_mars = {:f}'.format(e_mars))
    obliquity = 25.1893823  # degrees, Konopliv et al. 2016
    print('obliquity, degrees = {:f}'.format(obliquity))
    obliquity *= np.pi / 180

    n0 = 191.408  # degrees per year, Konopliv et al. 2011
    print('n0, degrees/yr = {:f}'.format(n0))
    n0 *= np.pi / 180 / 24 / 60 / 60 / 365.25

    C = (-1 / phi0dot) * (3/2) * (n0**2 / omega) * (1 - e_mars**2)**(-3/2) \
        * J2 * np.cos(obliquity)
    C *= mass_mars * r0_pot**2

    print('C (Konopliv et al. 2016) = {:e}'.format(C))
    print('C / (M R^2) (Konopliv et al. 2016) = {:e}'.format(C / mass_mars /
                                                             (r0_pot)**2))
    print('C / (M R_mpr^2) (Konopliv et al. 2016) = {:e}\n'
          .format(C / mass_mars / r_mars**2))

    II, AA, BB, CC, angles = InertiaTensor_from_C(C, potential, quiet=False,
                                                  normalize=True,
                                                  r_norm=r_mars)

    beta = 0.00032

    Acf = 0.0254088
    Bcf = 0.0254088
    Ccf = 0.0255172

    Ac = 0.0254036
    Bc = 0.0254113
    Cc = 0.0255199

    print('Free core nutation period (days), fluid planet = {:f}'
          .format(2 * np.pi / sigma_fcn(omega, AA, Acf, Bcf, Ccf, beta)
                  / 60 / 60 / 24))
    print('Free core nutation period (days), planet with lithosphere = {:f}'
          .format(2 * np.pi / sigma_fcn(omega, AA, Ac, Bc, Cc, beta)
                  / 60 / 60 / 24))


def sigma_fcn(omega, A, Ac, Bc, Cc, beta):
    '''
    Compute the free core nutation frequency.
    '''
    alpha = (2 * Cc - (Ac + Bc)) / (Ac + Bc)
    sigma = - omega * (A / (A - Ac)) * (alpha - beta)
    return sigma


# ==== EXECUTE SCRIPT ====


if __name__ == "__main__":
    main()
