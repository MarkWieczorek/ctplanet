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


CTPLANET
========

**ctplanet** provides several functions for working with the gravitational field of planetary crusts and hydrostatic density interfaces in the mantle and core. With this python package, you can:

* Compute the hydrostatic shape of density interfaces in an entirely hydrostatic planet.

* Compute the hydrostatic shape of density interfaces in a planet that has a non-hydrostatic lithosphere.

* Compute the relief along the crust-mantle interface from global gravity and topography data.

* Compute the inertia tensor of a planet when given the relief of all density interfaces.

* Read 1-D interior structure models in *deck* format.

In addition to these functions, example scripts are provided for computing crustal thickness models for the Moon and Mars, for computing the shapes of the fluid cores of the planets, and for computing free-core nutation periods.
