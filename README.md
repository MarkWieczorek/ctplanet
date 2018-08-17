# pyCrust
Create a crustal thickness map of a planet from gravity and topography.

## Description
The repo pyCrust provides several functions for generating crustal thickness maps of a planet from gravity and topography data, and the calculation of hydrostatic relief along density interfaces beneath the lithosphere.

### Modules
#### pyMoho.py
`pyMoho`                Calculate relief using a constant density crust and
                        mantle.

`pyMohoRho`             Calculate relief using a constant density mantle and a
                        variable density crust.
#### Hydrostatic.py
`HydrostaticShapeLith`  Calculate the relief of hydrostatic interfaces beneath
                        the lithosphere, along with the predicted gravity,
                        taking into account rotation and/or tides.

`HydrostaticShape`      Calculate the relief of hydrostatic interfaces and
                        predicted gravity of a rotating hydrostatic planet.

### Example scripts
`pyCrust_Moon`   A script that demonstrates how to calculate the thickenss
                 of the lunar crust using either a constant or variable density
                 crust. The latter can be used to reproduce the results
                 presented in *Wieczorek et al.* (2013).

`pyCrust_Mars`   A script that demonstrates how to calculate the thickenss
                 of the Martian crust using either a constant or variable
                 density crust. For the variable density crust, the density is
                 assumed to change discontinuously across the dichotomy
                 boundary.

`Core-Moon`      Calculate the hydrostatic relief of the lunar core accounting
                 for the non-hydrostatic potential that comes from the
                 lithosphere.

## How to install and run pyCrust

First, you will need to install [pyshtools](https://github.com/SHTOOLS/SHTOOLS), version 4.3 or greater. Normally, you should be able to do this with the
command

    pip install pyshtools

To install manually from the `develop` branch follow these steps

    git clone https://github.com/SHTOOLS/SHTOOLS.git
    cd SHTOOLS
    git checkout develop
    pip install -v -e .

Next install pyCrust in a different directory

    git clone
    https://github.com/MarkWieczorek/pycrust.git

Now run the script (preferably using python3)

    cd pycrust
    python pyCrust_Moon.py

## Reference

Wieczorek, M. A., G. A. Neumann, F. Nimmo, W. S. Kiefer, G. J. Taylor,
    H. J. Melosh, R. J. Phillips, S. C. Solomon, J. C. Andrews-Hanna,
    S. W. Asmar, A. S. Konopliv, F. G. Lemoine, D. E. Smith, M. M. Watkins,
    J. G. Williams, M. T. Zuber (2013), The crust of the Moon as seen by GRAIL,
    Science, 339, 671-675, doi:[10.1126/science.1231530](http://doi.org/10.1126/science.1231530).

