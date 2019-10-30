# pycrust
Create crustal thickness maps of planets from gravity and topography.

## Description
pycrust provides several functions and example scripts for generating crustal thickness maps of a planet from gravity and topography data, and the calculation of hydrostatic relief along density interfaces beneath the lithosphere.

### Methods
`pyMoho`                Calculate relief using a constant density crust and
                        mantle.

`pyMohoRho`             Calculate relief using a constant density mantle and a
                        variable density crust.

`HydrostaticShapeLith`  Calculate the relief of hydrostatic interfaces beneath
                        the lithosphere along with the predicted gravity,
                        taking into account rotation and/or tides using the
                        approach of *Wieczorek et al.* (2019).

`HydrostaticShape`      Calculate the relief of hydrostatic interfaces and
                        predicted gravity of a rotating hydrostatic planet
                        using the approach of *Wieczorek et al.* (2019).


`InertiaTensor_from_shape`  Calculate the inertia tensor given a radial density
                            profile and shape of each interface.

`InertiaTensor_from_C`      Calculate the inertia tensor given the polar moment
                            of inertia and the gravitational potential
                            coefficients.

`moi`                       Calculate the mean, normalized, moment of inertia
                            up to index n.


`ReadRefModel`   Read the reference interior model file.

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

`Mars_crust_thick_test`   Create a crustal thickness map of Mars from gravity
                          and topography and compare how results change if
                          hydrostatic interfaces are not taken into account.

`pyCrust_Mars_InSight`    Create a crustal thickness map of Mars from gravity
                          and topography, using the InSight crustal thickness
                          constraint.

`pyCrust_Mars_InSight_dichotomy`   Create a crustal thickness map of Mars from
                                   gravity and topography, using the InSight
                                   crustal thickness constraint and different
                                   densities across the dichotomy boundary.

`mars_fcn`       Compute the free core nutation period of Mars.

`mars_figs`      Create images related to Mars in *Wieczorek et al.* (2019).

`mars_j2`        Compute the contribution to the gravitational J2 of Mars from
                 hydrostatic interfaces beneath the lithosphere.

`Core-Moon`      Calculate the hydrostatic relief of the lunar core accounting
                 for the non-hydrostatic potential that comes from the
                 lithosphere.

`Earth_test`     Compute hydrostatic relief of Earth using PREM.

`ceres`          Calculate the hydrostatic shape of Ceres.

## How to install and run pyCrust

Download the pyCrust repository and install using pip
```bash
    git clone https://github.com/MarkWieczorek/pycrust.git
    pip install .
```

To execute a script
```bash
    cd examples
    python pyCrust_Moon.py
```

Depending on how your system is set up, it might be necessary to use explicitly `python3` and `pip3` instead of `python` and `pip` in the above commands.

## Reference

Wieczorek, M. A., G. A. Neumann, F. Nimmo, W. S. Kiefer, G. J. Taylor,
    H. J. Melosh, R. J. Phillips, S. C. Solomon, J. C. Andrews-Hanna,
    S. W. Asmar, A. S. Konopliv, F. G. Lemoine, D. E. Smith, M. M. Watkins,
    J. G. Williams, M. T. Zuber (2013), The crust of the Moon as seen by GRAIL,
    *Science*, 339, 671-675, doi:[10.1126/science.1231530](http://doi.org/10.1126/science.1231530).

Wieczorek, M. A.,  M. Beuthe, A. Rivoldini, and T. Van Hoolst (2019),
    Hydrostatic interfaces in bodies with nonhydrostatic lithospheres,
    *Journal of Geophysical Research: Planets*, 124, doi:[10.1029/2018JE005909](http://doi.org/10.1029/2018JE005909).
