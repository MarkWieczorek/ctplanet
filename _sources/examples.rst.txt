Examples
========

.. note::
    In order to access these example files, it will be necessary to either download or clone the entire ctplanet repo from `GitHub <https://github.com/MarkWieczorek/ctplanet>`_. The files will be located in the directory `examples`.

Moon
----

``Moon-Crust.py``
    A script that demonstrates how to calculate thethickenss of the lunar crust
    using either a constant or variable density crust. The latter can be used
    to reproduce the results presented in *Wieczorek et al.* (2013).

``Moon-Core.py``
    Calculate the hydrostatic relief of the lunar core accounting for the
    non-hydrostatic potential that comes from the lithosphere.

Mars
----

``Mars-Crust.py``
    A script that demonstrates how to calculate the thickenss of the Martian
    crust using either a constant or variable density crust. For the variable
    density crust, the density is assumed to change discontinuously across the
    dichotomy boundary.

``Mars-Crust-hydrostatic-tests.py``
    Create a crustal thickness map of Mars from gravity and topography and
    compare how results change if hydrostatic interfaces are not taken into
    account.

``Mars-Crust-InSight.py``
    Create a crustal thickness map of Mars from gravity and topography, using
    the InSight crustal thickness constraint.

``Mars-Crust-InSight-dichotomy.py``
    Create a crustal thickness map of Mars from gravity and topography, using
    the InSight crustal thickness constraint and different densities across
    the dichotomy boundary.

``Mars-fcn.py``
    Compute the free core nutation period of Mars.

``Mars-shape.py``
    Create images related to Mars in *Wieczorek et al.* (2019).

``Mars-j2.py``
    Compute the contribution to the gravitational J2 of Mars from hydrostatic
    interfaces beneath the lithosphere.

Earth
-----

``Earth-shape.py``
    Compute hydrostatic relief of Earth using PREM.

Ceres
-----

``Ceres-shape.py``
    Calculate the hydrostatic shape of Ceres.
