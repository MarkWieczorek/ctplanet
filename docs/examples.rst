Examples
========

.. note::
    In order to access the example files, it is necessary to download or clone the `ctplanet` repo from `GitHub <https://github.com/MarkWieczorek/ctplanet>`_. The files are located in the `examples` directory.

Moon
----

.. collapse:: <b>Moon-Crust.py</b>

    .. literalinclude:: ../examples/Moon-Crust.py

A script that demonstrates how to calculate the thickness of the lunar crust
using either a constant or variable density crust. The latter can be used
to reproduce the results presented in *Wieczorek et al.* (2013).

.. collapse:: <b>Moon-Core.py</b>

    .. literalinclude:: ../examples/Moon-Core.py

Calculate the hydrostatic relief of the lunar core accounting for the
non-hydrostatic potential that comes from the lithosphere.

Mars
----

.. collapse:: <b>Mars-Crust.py</b>

    .. literalinclude:: ../examples/Mars-Crust.py

A script that demonstrates how to calculate the thickness of the Martian
crust using either a constant or variable density crust. For the variable
density crust, the density is assumed to change discontinuously across the
dichotomy boundary.

.. collapse:: <b>Mars-Crust-hydrostatic-tests.py</b>

    .. literalinclude:: ../examples/Mars-Crust-hydrostatic-tests.py

Create a crustal thickness map of Mars from gravity and topography and
compare how results change if hydrostatic interfaces are not taken into
account.

.. collapse:: <b>Mars-Crust-InSight.py</b>

    .. literalinclude:: ../examples/Mars-Crust-InSight.py

Create a crustal thickness map of Mars from gravity and topography, using
the InSight crustal thickness constraint.

.. collapse:: <b>Mars-Crust-InSight-dichotomy.py</b>

    .. literalinclude:: ../examples/Mars-Crust-InSight-dichotomy.py

Create a crustal thickness map of Mars from gravity and topography, using
the InSight crustal thickness constraint and different densities across
the dichotomy boundary.

.. collapse:: <b>Mars-fcn.py</b>

    .. literalinclude:: ../examples/Mars-fcn.py

Compute the free core nutation period of Mars.

.. collapse:: <b>Mars-shape.py</b>

    .. literalinclude:: ../examples/Mars-shape.py

Create images related to Mars in *Wieczorek et al.* (2019).

.. collapse:: <b>Mars-j2.py</b>

    .. literalinclude:: ../examples/Mars-j2.py

Compute the contribution to the gravitational J2 of Mars from hydrostatic
interfaces beneath the lithosphere.

Earth
-----

.. collapse:: <b>Earth-shape.py</b>

    .. literalinclude:: ../examples/Earth-shape.py

Compute hydrostatic relief of Earth using PREM.

Ceres
-----

.. collapse:: <b>Ceres-shape.py</b>

    .. literalinclude:: ../examples/Ceres-shape.py

Calculate the hydrostatic shape of Ceres.
