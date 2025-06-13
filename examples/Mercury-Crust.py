#!/usr/bin/env python3
"""
Mercury-Crust

Create a crustal thickness map of Mercury from gravity and topography. 

The average crustal thickness is iterated in order to obtain a specified 
minimum crustal thickness.
"""

import os
import pyshtools as pysh

from ctplanet import pyMoho
from palettable import scientific as scm

# ==== MAIN FUNCTION ====


def main():

    directory = 'Mercury/'

    try:
        os.mkdir(directory)
    except:
        pass
        
    lmax_calc = 100
    lmax = 100
        
    potential = pysh.datasets.Mercury.JGMESS160A()
    
    print('Gravity file = {:s}'.format('JGMESS160A'))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}\n'.format(potential.gm))

    topo = pysh.datasets.Mercury.GTMES150(lmax=lmax)
    topo.r0 = topo.coeffs[0, 0, 0]

    print('Topography file = {:s}'.format('GTMES150'))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.r0 / 1.e3))
    
    filter = 1
    half = 40
    nmax = 7
    t0_sigma = 5.  # maximum difference between minimum crustal thickness
    
    t0 = 1.e3  # minimum crustal thickness
    print('Minimum crustal thickness (km) = {:f}'.format(t0/ 1.e3))

    rho_c = 2950.
    rho_m = 3200.
    print('Crustal density (kg/m3) = {:f}'.format(rho_c))
    print('Mantle density (kg/m3) = {:f}\n'.format(rho_m))
	
    tmin = 1.e9
    thickave = 35.e3  # initial guess of average crustal thickness

    while abs(tmin - t0) > t0_sigma:
        # iterate to fit assumed minimum crustal thickness

        moho = pyMoho(potential, topo, lmax, rho_c, rho_m,
                      thickave, filter_type=filter, half=half,
                      lmax_calc=lmax_calc, nmax=nmax, quiet=True)

        thick_grid = (topo.pad(lmax) - moho.pad(lmax)).expand(grid='DH2')
        print('Average crustal thickness (km) = {:f}'.format(thickave / 1.e3))
        tmin = thick_grid.min()
        tmax = thick_grid.max()
        print('Minimum thickness (km) = {:e}'.format(tmin / 1.e3))
        print('Maximum thickness (km) = {:e}'.format(tmax / 1.e3))
        thickave += t0 - tmin
                       
    # ---- Plot the crustal thickness grid and Moho spectrum ---
    print('\nPlotting the crustal thickness map and Moho spectrum ...')
          
    (thick_grid/1.e3).plot(show=False,
					cmap=scm.sequential.LaPaz_20.mpl_colormap,
					colorbar='right',
					cmap_limits=[0, 60],
					cb_label='Crustal thickness, km',
					cb_triangles='max',
					tick_interval=[60, 45],
					minor_tick_interval=[30, 15],
					fname=directory + 'Thick-Mercury.png')
					                 
    moho.plot_spectrum(show=False, fname=directory +'Moho-spectrum-Mercury.png')
    
	# ---- Save data to files ---
    print('Saving the crustal thickness grid to a file ...')
    (thick_grid/1.e3).to_file(directory + 'Thick-Mercury.dat')

    print('Saving the gridded data as a netcdf file for use with GMT ...')
    (thick_grid/1.e3).to_netcdf(directory + 'Thick-Mercury.nc')
    #(thick_grid).info()
    
    print("\nend program Mercury-Crust\n")
         
# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
