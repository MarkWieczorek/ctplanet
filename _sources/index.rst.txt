.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting Started

   installation.rst
   examples.rst
   references.rst

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference Documentation

   source/api.rst


.. figure:: Mars-crustal-thickness.png
    :width: 700px
    :align: center
    :alt: Crustal thickness of Mars

    Crustal thickness of Mars derived from gravity, topography, and InSight seismic constraints.

ctplanet
========

ctplanet provides several functions for working with the gravitational field of planetary crusts and hydrostatic density interfaces in the mantle and core. With this python module, you can:

* Compute the hydrostatic shape of density interfaces in an entirely hydrostatic planet.

* Compute the hydrostatic shape of density interfaces in a planet that has a non-hydrostatic lithosphere.

* Compute the relief along the crust-mantle interface from gravity and topography data.

* Compute the inertia tensor of a planet when given the relief of all density interfaces.

* Read 1-D interior structure models in 'deck' format.

In addition to these functions, example scripts are provided for computing crustal thickness models for the Moon and Mars, for computing the shapes of the fluid cores of the planets, and for computing free-core nutation periods.