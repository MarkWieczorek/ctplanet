# pyCrust
Create a crustal thickness map of a planet from gravity and topography.

## Description
`pyCrust` provides two functions for generating a crustal thickness map of a
planet from gravity and topography data. `pyMoho` assumes that the density of
the crust and mantle are both constant values, whereas `pyMohoRho` includes
the effects of lateral variations in crustal density.

The script `pyCrust_Moon` demonstrates how to use both functions, and can
be used to reproduce the results presented in *Wieczorek et al.* (2013).

This python program requires the use of [pyshtools](https://github.com/SHTOOLS/SHTOOLS), version 4.1 or greater.

## Reference
Wieczorek, M. A., G. A. Neumann, F. Nimmo, W. S. Kiefer, G. J. Taylor,
    H. J. Melosh, R. J. Phillips, S. C. Solomon, J. C. Andrews-Hanna,
    S. W. Asmar, A. S. Konopliv, F. G. Lemoine, D. E. Smith, M. M. Watkins,
    J. G. Williams, M. T. Zuber (2013), The crust of the Moon as seen by GRAIL,
    Science, 339, 671-675, doi:[10.1126/science.1231530](http://doi.org/10.1126/science.1231530).

