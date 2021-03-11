#!/usr/bin/env python3
"""
Venus-Crust-ave

Create a crustal thickness map of Venus from gravity and topography, assuming
an average crustal thickness. This script can be used to reproduce the results
presented in Wieczorek (2015).

Wieczorek, M. A. (2015). Gravity and topography of the terrestrial planets. 
	In T. Spohn, & G. Schubert (Eds.), Treatise on geophysics (second edition)
	Vol. 10, pp. 153–193. Oxford: Elsevier-Pergamon.
	https://doi.org/10.1016/B978-0-444-53802-4.00169-X
"""

import os
import pyshtools as pysh

from ctplanet import pyMoho
from palettable import scientific as scm

# ==== MAIN FUNCTION ====


def main():

    directory = 'Venus/average/'

    try:
        os.mkdir(directory)
    except:
        pass
    
    lmax_calc = 180
    lmax = lmax_calc * 2
    
    potential = pysh.datasets.Venus.MGNP180U()
    
    print('Gravity file = {:s}'.format('MGNP180U'))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}\n'.format(potential.gm))
    
    # read topo up to degree lmax
    topo = pysh.datasets.Venus.VenusTopo719(lmax=lmax)
    topo.r0 = topo.coeffs[0, 0, 0]

    print('Topography file = {:s}'.format('VenusTopo719'))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.r0 / 1.e3))
    
    filter = 1
    half = 70
    delta_max = 5.0
    nmax = 6
    thickave = 35.e3  # mean crustal thickness
    rho_c = 2900.0
    rho_m = 3330.0
    
    print('Average thickness of the crust (km) = {:e}'.format(thickave / 1.e3))
    print('Crustal density (kg/m3) = {:f}'.format(rho_c))
    print('Mantle density (kg/m3) = {:f}\n'.format(rho_m))
   
    moho = pyMoho(potential, topo, lmax, rho_c, rho_m, thickave,
                  filter_type=filter, half=half, lmax_calc=lmax_calc,
                   nmax=nmax, quiet=False, delta_max=delta_max)
  
    thick_grid = (topo.pad(lmax) - moho.pad(lmax)).expand(grid='DH2') / 1.e3
    
    # ---- Plot the crustal thickness grid and Moho spectrum ---
    print('\nPlotting the crustal thickness map and Moho spectrum ...')
          
    (thick_grid).plot(show=False,
					cmap=scm.sequential.LaPaz_20.mpl_colormap,
					colorbar='right',
					cmap_limits=[15, 75],
					cb_label='Crustal thickness, km',
					cb_triangles='both',
					tick_interval=[60, 45],
					minor_tick_interval=[30, 15],
					fname=directory + 'Thick-Venus-ave.png')
					                 
    moho.plot_spectrum(show=False, fname=directory +'Moho-spectrum-Venus-ave.png')
    
	# ---- Save data to files ---
    print('Saving the crustal thickness grid to a file ...')
    (thick_grid).to_file(directory + 'Thick-Venus-ave.dat')

    print('Saving the gridded data as a netcdf file for use with GMT ...')
    (thick_grid).to_netcdf(directory + 'Thick-Venus-ave.nc')
    #(thick_grid).info()
    
    print("\nend program Venus-Crust-ave\n")
         
# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
