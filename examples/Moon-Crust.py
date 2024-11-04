#!/usr/bin/env python3
"""
Moon-Crust

Create a crustal thickness map of the Moon from gravity and topography.

This script generates two different crustal thickness maps. The first assumes
that the density of both the crust and mantle are constant, whereas the second
includes the effect of lateral variations in crustal density. This script can
be used to reproduce the results presented in Wieczorek et al. (2013).

Wieczorek, M. A., G. A. Neumann, F. Nimmo, W. S. Kiefer, G. J. Taylor,
    H. J. Melosh, R. J. Phillips, S. C. Solomon, J. C. Andrews-Hanna,
    S. W. Asmar, A. S. Konopliv, F. G. Lemoine, D. E. Smith, M. M. Watkins,
    J. G. Williams, M. T. Zuber (2013), The crust of the Moon as seen by GRAIL,
    Science, 339, 671-675, doi:10.1126/science.1231530, 2013.
"""
import os
import pyshtools as pysh

from ctplanet import pyMoho
from ctplanet import pyMohoRho

# ==== MAIN FUNCTION ====


def main():

    lmax = 900  # determines spatial resolution of grids
    lmax_calc = 600  # maximum degree to use in calculations

    densityfile = 'Data/density_no_mare_n3000_f3050_719.sh'

    pot = pysh.datasets.Moon.GRGM900C()

    try:
        os.mkdir('figs')
    except Exception:
        pass

    print('Gravity file = {:s}'.format('GRGM900C'))
    print('Lmax of potential coefficients = {:d}'.format(pot.lmax))
    print('Reference radius (m) = {:e}'.format(pot.r0))
    print('GM = {:e}\n'.format(pot.gm))

    topo = pysh.datasets.Moon.LDEM_shape_pa(lmax=lmax)
    topo.r0 = topo.coeffs[0, 0, 0]

    print('Topography file = {:s}'.format('LDEM_shape_pa'))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (m) = {:e}\n'.format(topo.r0))

    density = pysh.SHCoeffs.from_file(densityfile, lmax=719)
    rho_c0 = density.coeffs[0, 0, 0]  # average grain density

    print('Average grain density of crust (kg/m3) = {:e}'.format(rho_c0))
    print('Lmax of density coefficients = {:d}\n'.format(density.lmax))

    a_12_14_lat = -3.3450
    a_12_14_long = -20.450

    # Parameters corresponding to model 1 of Wieczorek et al. (2013),
    # with the exception that the average thickness has been increased by 1 km.
    thickave = 35.e3
    porosity = 0.12
    rho_m = 3220.0
    filter = 1
    half = 80

    print('Average thickness of the crust (km) = {:e}'.format(thickave / 1.e3))
    print('Porosity (%)= {:e}'.format(porosity * 100))
    print('Mantle density (kg/m3)= {:e}'.format(rho_m))

    # Constant density model
    print('\n=== Constant density crust ===')
    rho_c = rho_c0 * (1. - porosity)  # assumed constant density
    print('Bulk density of the crust(kg/m3)= {:e}'.format(rho_c * 100))

    moho = pyMoho(pot, topo, lmax, rho_c, rho_m, thickave,
                  filter_type=filter, half=half, lmax_calc=lmax_calc,
                  quiet=False, delta_max=50.)

    thick_grid = (topo.pad(lmax) - moho.pad(lmax)).expand(grid='DH2') / 1.e3
    thick_grid.plot(show=False, colorbar='bottom',
                    cb_label='Crustal thickness, km',
                    fname='figs/Thick-Moon-1.png')
    moho.plot_spectrum(show=False, fname='figs/Moho-spectrum-Moon-1.png')

    print('Crustal thickness at Apollo 12/14 landing sites (km) = {:e}'
          .format((topo.pad(lmax) - moho.pad(lmax))
                  .expand(lat=a_12_14_lat, lon=a_12_14_long) / 1.e3))

    # Model with variable density
    print('\n=== Variable density crust ===')

    moho = pyMohoRho(pot, topo, density, porosity, lmax, rho_m,
                     thickave, filter_type=filter, half=half,
                     lmax_calc=lmax_calc, quiet=False, delta_max=50.)

    thick_grid = (topo-moho.pad(topo.lmax)).expand(grid='DH2') / 1.e3
    thick_grid.plot(show=False, colorbar='bottom',
                    cb_label='Crustal thickness, km',
                    fname='figs/Thick-Moon-2.png')
    moho.plot_spectrum(show=False, fname='figs/Moho-spectrum-Moon-2.png')

    print('Crustal thickness at Apollo 12/14 landing sites (km) = {:e}'.format(
        (topo-moho.pad(topo.lmax)).expand(lat=a_12_14_lat,
                                          lon=a_12_14_long) / 1.e3))


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
