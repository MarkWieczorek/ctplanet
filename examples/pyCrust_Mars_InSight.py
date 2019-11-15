#!/usr/bin/env python3
"""
pyCrust_Mars_InSight

Create a crustal thickness map of Mars from gravity and topography, using
the InSight crustal thickness constraint.

The script assumes that the density of both the crust and mantle are constant,
and that porosity in present in a constant thickness surficial layer.
The average crustal thickness is iterated in order to fit the specified
thickness at the InSight landing site.
"""
import os
import matplotlib.pyplot as plt

import pyshtools

from pycrust import pyMoho
from pycrust import HydrostaticShapeLith
from pycrust import ReadRefModel

# ==== MAIN FUNCTION ====


def main():

    gravfile = 'Data/gmm3_120_sha.tab'
    topofile = 'Data/MarsTopo719.shape'

    model_name = ['DWThot', 'DWThotCrust1', 'DWThotCrust1r', 'EH45Tcold',
                  'EH45TcoldCrust1', 'EH45TcoldCrust1r', 'EH45ThotCrust2',
                  'EH45ThotCrust2r', 'LFAK', 'SANAK', 'TAYAK', 'DWAK',
                  'ZG_DW']
    spec = 'Data/Mars-reference-interior-models/Smrekar/'
    interior_file = [spec + name + '.deck' for name in model_name]

    lmax_calc = 90
    lmax = lmax_calc * 4

    potential = pyshtools.SHGravCoeffs.from_file(gravfile, header_units='km')

    directory = 'InSight/constant-porosity/'

    try:
        os.mkdir(directory)
    except:
        pass

    print('Gravity file = {:s}'.format(gravfile))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}\n'.format(potential.gm))

    topo = pyshtools.SHCoeffs.from_file(topofile, lmax=lmax)
    topo.r0 = topo.coeffs[0, 0, 0]

    print('Topography file = {:s}'.format(topofile))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.r0 / 1.e3))

    filter = 1
    half = 50
    nmax = 7
    lmax_hydro = 15
    t_sigma = 5.  # maximum difference between crustal thickness iterations
    omega = pyshtools.constant.omega_mars.value

    d_lith = 150.e3

    # Get input model and InSight crustal thickness

    lat_insight = 4.502384
    lon_insight = 135.623447
    print('Input crustal thickness at InSight landing site (km) >')
    t_insight_ref = float(input())
    t_insight_ref *= 1.e3

    print('Input thickness of porous layer (km) >')
    t_porosity = float(input())
    t_porosity *= 1.e3

    porosity = 0.
    if t_porosity != 0.:
        print('Input porosity of porous layer in percent >')
        porosity = float(input())
        porosity /= 100.

    rho_c_min = 2550.
    rho_c_max = 3200.
    rho_c_int = 50.

    f_summary = open(directory + 'summary.txt', 'w')
    f_summary.write('Model    rho_c    rho_mantle    t_insight    t_ave    '
                    + ' t_min    t_max\n')

    # Compute contribution to the gravitational potential resulting
    # from porosity. The density contrast at the surface of this layer
    # is: - porosity * rho_crust and at the base it is porosity * rho_crust.
    # Here we assume that rho_c is 1 and then multiply this later by the
    # correct value.
    if porosity != 0.:
        pot_porosity = pyshtools.SHGravCoeffs.from_shape(topo, -porosity,
                                                         potential.gm,
                                                         nmax=nmax,
                                                         lmax_grid=lmax,
                                                         lmax=lmax_calc)
        pot_porosity = pot_porosity.change_ref(r0=potential.r0)

        pot_lower = pyshtools.SHGravCoeffs.from_shape(topo - t_porosity,
                                                      porosity,
                                                      potential.gm,
                                                      nmax=nmax,
                                                      lmax_grid=lmax,
                                                      lmax=lmax_calc)
        pot_porosity += pot_lower.change_ref(r0=potential.r0)

    for model in range(len(model_name)):
        print('Working on model : {:s}'.format(model_name[model]))

        # --- read 1D reference interior model ---
        radius, rho, i_crust, i_core, i_lith = ReadRefModel(
            interior_file[model], depth=d_lith, quiet=False)

        rho_mantle = rho[i_crust-1]
        rho_core = rho[i_core-1]
        # n = len(radius) - 1
        # r0_model = radius[n]

        for i in range(int((rho_c_max - rho_c_min) / rho_c_int + 1)):
            rho_c = rho_c_min + i * rho_c_int
            print('rho_crust (kg/m3) = {:f}'.format(rho_c))

            ident = model_name[model] + '-' + str(int(t_insight_ref / 1.e3)) \
                + '-' + str(int(rho_c))
            if porosity != 0.:
                ident += '-' + str(int(t_porosity / 1.e3)) + '-' \
                    + str(int(porosity * 100))

            # Compute gravity contribution from hydrostatic density interfaces
            thickave = 44.e3  # initial guess of average crustal thickness
            r_sigma = topo.r0 - thickave

            hlm, clm_hydro, mass_model = \
                HydrostaticShapeLith(radius, rho, i_lith, potential, topo,
                                     rho_c, r_sigma, omega, lmax_hydro)

            print('Total mass of model (kg) = {:e}'.format(mass_model))
            print('% of J2 arising from beneath lithosphere = {:f}'
                  .format(clm_hydro.coeffs[0, 2, 0] /
                          potential.coeffs[0, 2, 0] * 100.))

            pot_lith = potential.copy()
            pot_lith.coeffs[:, :lmax_hydro+1, :lmax_hydro+1] -= \
                clm_hydro.coeffs[:, :lmax_hydro+1, :lmax_hydro+1]

            if porosity != 0.:
                pot_lith.coeffs[:, :lmax_calc+1, :lmax_calc+1] -= \
                    rho_c * pot_porosity.coeffs[:, :lmax_calc+1, :lmax_calc+1]

            t_insight = 1.e9
            thickave = 44.e3    # initial guess of average crustal thickness

            pot_temp = pot_lith.copy()
            pot_porosity_correction = pyshtools.SHGravCoeffs.from_zeros(
                lmax_calc, potential.gm, potential.r0)

            while abs(t_insight_ref - t_insight) > t_sigma:
                # iterate to fit assumed minimum crustal thickness

                if porosity != 0.:
                    pot_temp = pot_lith.pad(lmax=lmax_calc) \
                        - pot_porosity_correction.pad(lmax=lmax_calc)

                moho = pyMoho(pot_temp, topo, lmax, rho_c, rho_mantle,
                              thickave, filter_type=filter, half=half,
                              lmax_calc=lmax_calc, nmax=nmax,
                              quiet=True)

                if porosity != 0.:
                    topo_grid = topo.pad(lmax).expand(grid='DH2')
                    moho_grid = moho.pad(lmax).expand(grid='DH2')
                    thick_grid = topo_grid - moho_grid
                else:
                    thick_grid = (topo.pad(lmax)
                                  - moho.pad(lmax)).expand(grid='DH2')
                t_insight = (topo.pad(lmax) -
                             moho.pad(lmax)).expand(lat=lat_insight,
                                                    lon=lon_insight)

                tmin = thick_grid.min()
                tmax = thick_grid.max()
                thickave += t_insight_ref - t_insight

                # Compute correction resulting from porosity in the mantle
                if porosity != 0.:
                    lower_grid = topo_grid - t_porosity
                    upper_grid = lower_grid.copy()
                    upper_grid.data[moho_grid.data > lower_grid.data] = \
                        moho_grid.data[moho_grid.data > lower_grid.data]

                    pot_porosity_correction = \
                        pyshtools.SHGravCoeffs.from_shape(
                            upper_grid, -porosity * (rho_mantle - rho_c),
                            potential.gm, nmax=nmax, lmax_grid=lmax,
                            lmax=lmax_calc)
                    pot_porosity_correction = \
                        pot_porosity_correction.change_ref(r0=potential.r0)

                    pot_lower = pyshtools.SHGravCoeffs.from_shape(
                        lower_grid, porosity * (rho_mantle - rho_c),
                        potential.gm, nmax=nmax, lmax_grid=lmax,
                        lmax=lmax_calc)
                    pot_porosity_correction += \
                        pot_lower.change_ref(r0=potential.r0)

            print('Average crustal thickness (km) = {:f}'
                  .format(thickave / 1.e3))
            print('Crustal thickness at InSight landing sites (km) = {:f}'
                  .format(t_insight / 1.e3))
            print('Minimum thickness (km) = {:e}'.format(tmin / 1.e3))
            print('Maximum thickness (km) = {:e}'.format(tmax / 1.e3))

            if tmin < 0:
                break

            moho.pad(lmax_calc).to_file(directory + 'Moho-Mars-'
                                        + ident + '.sh')
            fig, ax = (thick_grid/1.e3).plot(
                show=False, colorbar=True,
                cb_label='Crustal thickness (km)',
                cb_orientation='horizontal',
                fname=directory + 'Thick-Mars-' + ident + '.png')
            # fig2, ax2 = moho.plot_spectrum(
            #    show=False,
            #    fname=directory + 'Moho-spectrum-Mars-' + ident + '.png',
            #    legend=ident)
            f_summary.write('{:s}    {:f}    {:f}    {:f}    {:f}    '
                            '{:f}    {:f}\n'
                            .format(model_name[model], rho_c, rho_mantle,
                                    t_insight/1.e3, thickave / 1.e3,
                                    tmin / 1.e3, tmax / 1.e3))
            plt.close(fig)
            # plt.close(fig2)

    f_summary.close()


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
