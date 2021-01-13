Examples
========

The following example scripts can be found in the directory `examples`:

``pyCrust_Moon.py``

    A script that demonstrates how to calculate thethickenss of the lunar crust
    using either a constant or variable density crust. The latter can be used
    to reproduce the results presented in *Wieczorek et al.* (2013).

``pyCrust_Mars.py``

    A script that demonstrates how to calculate the thickenss of the Martian
    crust using either a constant or variable density crust. For the variable
    density crust, the density is assumed to change discontinuously across the
    dichotomy boundary.

``Mars_crust_thick_test.py``

    Create a crustal thickness map of Mars from gravity and topography and
    compare how results change if hydrostatic interfaces are not taken into
    account.

``pyCrust_Mars_InSight.py``

    Create a crustal thickness map of Mars from gravity and topography, using
    the InSight crustal thickness constraint.

``pyCrust_Mars_InSight_dichotomy.py``

    Create a crustal thickness map of Mars from gravity and topography, using
    the InSight crustal thickness constraint and different densities across
    the dichotomy boundary.

``mars_fcn.py``

    Compute the free core nutation period of Mars.

``mars_figs.py``

    Create images related to Mars in *Wieczorek et al.* (2019).

``mars_j2.py``

    Compute the contribution to the gravitational J2 of Mars from hydrostatic
    interfaces beneath the lithosphere.

``Core-Moon.py``

    Calculate the hydrostatic relief of the lunar core accounting for the
    non-hydrostatic potential that comes from the lithosphere.

``Earth_test.py``

    Compute hydrostatic relief of Earth using PREM.

``ceres.py``

    Calculate the hydrostatic shape of Ceres.
