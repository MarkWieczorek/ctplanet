#!/usr/bin/env python3
"""
Mars-Crust-InSight-dichotomy

Create a crustal thickness map of Mars from gravity and topography, using
the InSight crustal thickness constraint.

The script assumes that the mantle density is constant and that the density of
the crust differs on either side of the dichotomy boundary. The gravitational
contribution of the polar caps are explicitly accounted for, and the average
crustal thickness is iterated in order to fit the specified thickness at the
InSight landing site.
"""
import os
import matplotlib.pyplot as plt
import pyshtools as pysh

from ctplanet import pyMohoRho
from ctplanet import HydrostaticShapeLith
from ctplanet import ReadRefModel

# ==== MAIN FUNCTION ====


def main():

    densityfile = 'Data/dichotomy_359.sh'
    northpolarcap = 'Data/Mars_NorthPolarCapThickness719.sh.gz'
    southpolarcap = 'Data/Mars_SouthPolarCapThickness719.sh.gz'

    model_name1 = ['DWThot', 'DWThotCrust1', 'DWThotCrust1r', 'EH45Tcold',
                   'EH45TcoldCrust1', 'EH45TcoldCrust1r', 'EH45ThotCrust2',
                   'EH45ThotCrust2r', 'LFAK', 'SANAK', 'TAYAK', 'DWAK',
                   'ZG_DW', 'YOTHotRc1760kmDc40km', 'YOTHotRc1810kmDc40km']
    spec1 = 'Data/Mars-reference-interior-models/Smrekar/'
    model_name2 = ['Khan2022']
    spec2 = 'Data/Mars-reference-interior-models/'

    interior_file = [spec1 + name + '.deck' for name in model_name1]
    interior_file += [spec2 + name + '.deck' for name in model_name2]
    model_name = model_name1 + model_name2

    lmax_calc = 90
    lmax = lmax_calc * 4

    potential = pysh.datasets.Mars.GMM3()

    directory = 'InSight/dichotomy/'

    try:
        os.mkdir(directory)
    except Exception:
        pass

    print('Gravity file = {:s}'.format('GMM3'))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}\n'.format(potential.gm))

    topo = pysh.datasets.Mars.MOLA_shape(lmax=lmax)
    topo_grid = topo.pad(lmax).expand(grid='DH2')

    print('Topography file = {:s}'.format('MOLA_shape'))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.coeffs[0, 0, 0] / 1.e3))

    print('Crustal density file = {:s}'.format(densityfile))
    dichotomy = pysh.SHCoeffs.from_file(densityfile).pad(lmax=lmax)

    print('Lmax of density coefficients = {:d}\n'.format(dichotomy.lmax))

    # read polar cap thicknesses up to lmax
    npc = pysh.SHCoeffs.from_file(northpolarcap, lmax=lmax)
    spc = pysh.SHCoeffs.from_file(southpolarcap, lmax=lmax)
    rho_npc = 1250.
    rho_spc = 1300.

    print('North polar cap thickness file = {:s}'.format(northpolarcap))
    print('Lmax of North polar cap coefficients = {:d}'.format(npc.lmax))
    print('Assumed density of the North polar cap (kg m-3) = {:e}'
          .format(rho_npc))
    print('South polar cap thickness file = {:s}'.format(southpolarcap))
    print('Lmax of South polar cap coefficients = {:d}'.format(spc.lmax))
    print('Assumed density of the South polar cap (kg m-3) = {:e}'
          .format(rho_npc))

    # Topography excluding the polar caps
    topo_no_pc = topo - npc - spc

    # compute gravitational contribution of the polar caps
    npc_correction = pysh.SHGravCoeffs.from_shape(shape=topo_no_pc+npc,
                                                  rho=1.0,
                                                  gm=potential.gm,
                                                  nmax=8,
                                                  lmax=lmax_calc,
                                                  lmax_grid=lmax,
                                                  lmax_calc=lmax)
    spc_correction = pysh.SHGravCoeffs.from_shape(shape=topo_no_pc + spc,
                                                  rho=1.0,
                                                  gm=potential.gm,
                                                  nmax=8,
                                                  lmax=lmax_calc,
                                                  lmax_grid=lmax,
                                                  lmax_calc=lmax)
    base_correction = pysh.SHGravCoeffs.from_shape(shape=topo_no_pc,
                                                   rho=1.0,
                                                   gm=potential.gm,
                                                   nmax=8,
                                                   lmax=lmax_calc,
                                                   lmax_grid=lmax,
                                                   lmax_calc=lmax)
    pc_correction = rho_npc * npc_correction.change_ref(r0=potential.r0)
    pc_correction += rho_spc * spc_correction.change_ref(r0=potential.r0)
    pc_correction -= (rho_spc + rho_npc) * \
        base_correction.change_ref(r0=potential.r0)

    filter = 1
    half = 50
    nmax = 7
    lmax_hydro = 15
    t_sigma = 5.  # maximum difference between crustal thickness iterations
    omega = pysh.constants.Mars.angular_velocity.value

    d_lith = 150.e3

    # Get input model and InSight crustal thickness

    lat_insight = 4.502384
    lon_insight = 135.623447
    print('Input crustal thickness at InSight landing site (km) >')
    t_insight_ref = float(input())
    t_insight_ref *= 1.e3

    rho_c_min = 2550.
    rho_c_max = 3200.
    rho_c_int = 50.

    ident = str(int(t_insight_ref / 1.e3))
    f_summary = open(directory + 'summary-' + ident + '.txt', 'w')
    f_summary.write('Model    rho_south    rho_north    rho_mantle    '
                    't_insight    t_ave    t_min    t_max\n')

    for model in range(len(model_name)):
        print('Working on model : {:s}'.format(model_name[model]))

        # --- read 1D reference interior model ---
        radius, rho, i_crust, i_core, i_lith = ReadRefModel(
            interior_file[model], depth=d_lith, quiet=False)

        rho_mantle = rho[i_crust-1]
        # rho_core = rho[i_core-1]
        # n = len(radius) - 1
        # r0_model = radius[n]

        for i in range(int((rho_c_max - rho_c_min) / rho_c_int + 1)):
            rho_south = rho_c_min + i * rho_c_int

            for j in range(int((rho_c_max - rho_c_min) / rho_c_int + 1)):
                rho_north = rho_c_min + j * rho_c_int

                if rho_north < rho_south:
                    continue

                print('rho_south = {:f}\n'.format(rho_south) +
                      'rho_north = {:f}'.format(rho_north))

                density = dichotomy.copy()
                density = density * (rho_north - rho_south)
                density.coeffs[0, 0, 0] += rho_south
                crustal_porosity = 0.0

                ident = model_name[model] + '-' \
                    + str(int(t_insight_ref / 1.e3)) \
                    + '-' + str(int(rho_south)) + '-' + str(int(rho_north))

                # Compute gravity contribution from hydrostatic density
                # interfaces
                thickave = 44.e3  # initial guess of average crustal thickness
                r_sigma = topo.coeffs[0, 0, 0] - thickave

                hlm, clm_hydro, mass_model = \
                    HydrostaticShapeLith(radius, rho, i_lith, potential, topo,
                                         (rho_south + rho_north) / 2, r_sigma,
                                         omega, lmax_hydro)

                print('Total mass of model (kg) = {:e}'.format(mass_model))
                print('% of J2 arising from beneath lithosphere = {:f}'
                      .format(clm_hydro.coeffs[0, 2, 0] /
                              potential.coeffs[0, 2, 0] * 100.))

                pot_lith = potential.copy()
                pot_lith.coeffs[:, :lmax_hydro+1, :lmax_hydro+1] -= \
                    clm_hydro.coeffs[:, :lmax_hydro+1, :lmax_hydro+1]

                t_insight = 1.e9
                thickave = 44.e3    # initial guess of average thickness

                while abs(t_insight_ref - t_insight) > t_sigma:
                    # iterate to fit assumed minimum crustal thickness

                    moho = pyMohoRho(pot_lith, topo_no_pc, density,
                                     crustal_porosity, lmax,
                                     rho_mantle, thickave,
                                     filter_type=filter,
                                     half=half, lmax_calc=lmax_calc,
                                     quiet=True,
                                     correction=pc_correction,
                                     nmax=nmax)
                    moho_grid = moho.pad(lmax).expand(grid='DH2')

                    thick_grid = topo_grid - moho_grid
                    t_insight = (topo.pad(lmax) -
                                 moho.pad(lmax)).expand(lat=lat_insight,
                                                        lon=lon_insight)

                    tmin = thick_grid.min()
                    tmax = thick_grid.max()
                    thickave += t_insight_ref - t_insight

                print('Average crustal thickness (km) = {:f}'
                      .format(thickave / 1.e3))
                print('Crustal thickness at InSight landing sites (km) = {:f}'
                      .format(t_insight / 1.e3))
                print('Minimum thickness (km) = {:e}'.format(tmin / 1.e3))
                print('Maximum thickness (km) = {:e}'.format(tmax / 1.e3))

                if tmin < 0:
                    break

                moho.pad(lmax_calc).to_file(
                    directory + 'Moho-Mars-' + ident + '.sh')
                fig, ax = (thick_grid/1.e3).plot(
                    show=False,
                    colorbar='bottom',
                    cb_label='Crustal thickness (km)',
                    fname=directory + 'Thick-Mars-' + ident + '.png')
                # fig2, ax2 = moho.plot_spectrum(
                #    show=False,
                #    fname=directory + 'Moho-spectrum-Mars-'
                #    + ident + '.png',
                #    legend=ident)
                f_summary.write('{:s}    {:f}    {:f}    {:f}    {:f}    '
                                '{:f}    {:f}    {:f}\n'
                                .format(model_name[model], rho_south,
                                        rho_north, rho_mantle,
                                        t_insight / 1.e3,
                                        thickave / 1.e3, tmin / 1.e3,
                                        tmax / 1.e3))
                plt.close(fig)
                # plt.close(fig2)

    f_summary.close()


# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
