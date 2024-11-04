#!/usr/bin/env python3
"""
Mars-Crust-InSight

Create a crustal thickness map of Mars from gravity and topography, using
the InSight crustal thickness constraint.

The script assumes that the density of both the crust and mantle are constant,
and that porosity in present in a constant thickness surficial layer (excluding
the polar caps). The gravitational contribution of the polar caps are
explicitly accounted for, and the average crustal thickness is iterated in
order to fit the specified thickness at the InSight landing site.
"""
import os
import matplotlib.pyplot as plt
import pyshtools as pysh

from ctplanet import pyMoho
from ctplanet import HydrostaticShapeLith
from ctplanet import ReadRefModel

# ==== MAIN FUNCTION ====


def main():

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

    print('Gravity file = {:s}'.format('GMM3'))
    print('Lmax of potential coefficients = {:d}'.format(potential.lmax))
    print('Reference radius (km) = {:f}'.format(potential.r0 / 1.e3))
    print('GM = {:e}\n'.format(potential.gm))

    # read topo up to degree lmax
    topo = pysh.datasets.Mars.MOLA_shape(lmax=lmax)
    topo_grid = topo.pad(lmax).expand(grid='DH2')

    print('Topography file = {:s}'.format('MOLA_shape'))
    print('Lmax of topography coefficients = {:d}'.format(topo.lmax))
    print('Reference radius (km) = {:f}\n'.format(topo.coeffs[0, 0, 0] / 1.e3))

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
    topo_no_pc_grid = topo_no_pc.expand(grid='DH2')

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
    ident = str(int(t_insight_ref / 1.e3))

    print('Input thickness of porous layer (km) >')
    t_porosity = float(input())
    t_porosity *= 1.e3

    porosity = 0.
    if t_porosity != 0.:
        directory = 'InSight/constant-porosity/'
        print('Input porosity of porous layer in percent >')
        porosity = float(input())
        porosity /= 100.
    else:
        directory = 'InSight/constant/'

    if porosity != 0.:
        ident += '-' + str(int(t_porosity / 1.e3)) + '-' \
            + str(int(porosity * 100))

    try:
        os.mkdir(directory)
    except Exception:
        pass

    rho_c_min = 2550.
    rho_c_max = 3300.
    rho_c_int = 50.

    f_summary = open(directory + 'summary-' + ident + '.txt', 'w')
    f_summary.write('Model    rho_c    rho_mantle    t_insight    t_ave    '
                    + ' t_min    t_max\n')

    # Compute contribution to the gravitational potential resulting
    # from porosity. The density contrast at the surface of this layer
    # is: - porosity * rho_crust and at the base it is porosity * rho_crust.
    # Here we assume that rho_c is 1 and then multiply this later by the
    # correct value. Note that no porosity is added to the polar caps.
    if porosity != 0.:
        pot_porosity = pysh.SHGravCoeffs.from_shape(topo_no_pc, -porosity,
                                                    potential.gm,
                                                    nmax=nmax,
                                                    lmax=lmax_calc,
                                                    lmax_grid=lmax,
                                                    lmax_calc=lmax
                                                    )
        pot_porosity = pot_porosity.change_ref(r0=potential.r0)

        pot_lower = pysh.SHGravCoeffs.from_shape(topo_no_pc - t_porosity,
                                                 porosity,
                                                 potential.gm,
                                                 nmax=nmax,
                                                 lmax=lmax_calc,
                                                 lmax_grid=lmax,
                                                 lmax_calc=lmax
                                                 )
        pot_porosity += pot_lower.change_ref(r0=potential.r0)

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
            rho_c = rho_c_min + i * rho_c_int
            print('rho_crust (kg/m3) = {:f}'.format(rho_c))

            ident = model_name[model] + '-' + str(int(t_insight_ref / 1.e3)) \
                + '-' + str(int(rho_c))
            if porosity != 0.:
                ident += '-' + str(int(t_porosity / 1.e3)) + '-' \
                    + str(int(porosity * 100))

            # Compute gravity contribution from hydrostatic density interfaces
            thickave = 44.e3  # initial guess of average crustal thickness
            r_sigma = topo.coeffs[0, 0, 0] - thickave

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

            pot_temp = pot_lith.copy()
            pot_porosity_correction = pysh.SHGravCoeffs.from_zeros(
                lmax_calc, potential.gm, potential.r0)

            # iterate to fit assumed InSight crustal thickness
            while abs(t_insight_ref - t_insight) > t_sigma:
                if porosity != 0.:
                    pot_temp = pot_lith.pad(lmax=lmax_calc) \
                        - pot_porosity_correction.pad(lmax=lmax_calc)

                moho = pyMoho(pot_temp, topo_no_pc, lmax, rho_c, rho_mantle,
                              thickave, filter_type=filter, half=half,
                              lmax_calc=lmax_calc, nmax=nmax,
                              correction=pc_correction, quiet=True)
                moho_grid = moho.pad(lmax).expand(grid='DH2')

                thick_grid = topo_grid - moho_grid
                t_insight = (topo.pad(lmax) -
                             moho.pad(lmax)).expand(lat=lat_insight,
                                                    lon=lon_insight)

                tmin = thick_grid.min()
                tmax = thick_grid.max()
                thickave += t_insight_ref - t_insight

                # Compute correction resulting from porosity in the mantle
                if porosity != 0.:
                    lower_grid = topo_no_pc_grid - t_porosity
                    upper_grid = lower_grid.copy()
                    upper_grid.data[moho_grid.data > lower_grid.data] = \
                        moho_grid.data[moho_grid.data > lower_grid.data]

                    pot_porosity_correction = \
                        pysh.SHGravCoeffs.from_shape(
                            upper_grid, -porosity * (rho_mantle - rho_c),
                            potential.gm, nmax=nmax, lmax_grid=lmax,
                            lmax=lmax_calc)
                    pot_porosity_correction = \
                        pot_porosity_correction.change_ref(r0=potential.r0)

                    pot_lower = pysh.SHGravCoeffs.from_shape(
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
                show=False, colorbar='bottom',
                cb_label='Crustal thickness (km)',
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
