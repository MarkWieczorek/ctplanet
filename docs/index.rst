pyCrust
=======

pyCrust provides several functions and example scripts for generating
crustal thickness maps of a planet from gravity and topography data. 

How to install and run pyCrust
------------------------------

Download the pyCrust repository and install using pip

.. code:: bash

    git clone https://github.com/MarkWieczorek/pycrust.git
    pip install .

To execute a script

.. code:: bash

    cd examples
    python pyCrust_Moon.py

Depending on how your system is set up, it might be necessary to use
explicitly ``python3`` and ``pip3`` instead of ``python`` and ``pip`` in
the above commands.

References
----------

Wieczorek, M. A., G. A. Neumann, F. Nimmo, W. S. Kiefer, G. J. Taylor,
H. J. Melosh, R. J. Phillips, S. C. Solomon, J. C. Andrews-Hanna, S. W.
Asmar, A. S. Konopliv, F. G. Lemoine, D. E. Smith, M. M. Watkins, J. G.
Williams, M. T. Zuber (2013), The crust of the Moon as seen by GRAIL,
*Science*, 339, 671-675,
doi:\ `10.1126/science.1231530 <http://doi.org/10.1126/science.1231530>`__.

Wieczorek, M. A., M. Beuthe, A. Rivoldini, and T. Van Hoolst (2019),
Hydrostatic interfaces in bodies with nonhydrostatic lithospheres,
*Journal of Geophysical Research: Planets*, 124,
doi:\ `10.1029/2018JE005909 <http://doi.org/10.1029/2018JE005909>`__.


.. toctree::
   :maxdepth: 2
   :caption: Documentation:

   routines.rst
   examples.rst
