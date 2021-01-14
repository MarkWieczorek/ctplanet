'''
Read the reference interior model file.
'''
import numpy as np


# ==== ReadRefModel ====

def ReadRefModel(filename, depth=None, quiet=True):
    '''
    Read a reference interior model file in deck format and convert to the form
    required by HydrostaticShapeLith.

    Returns
    -------
    radius : ndarray, size (n+1)
        A vector of radii, where radius[0] is the center of the planet and
        radius[n] is the surface.
    rho : ndarray, size (n+1)
        A vector of densities of the layers, where rho[i] is the density
        between radius[i] and radius[i+1]. The density above the surface,
        rho[n], is set to zero. The density of the base of the crust is
        rho[i_crust], the density of the upper mantle is rho[i_crust-1], and
        the density of the upper core is rho[i_core-1].
    i_crust : integer
        Index of the radius vector corresponding to the base of the crust.
    i_core : integer
        Index of the radius vector corresponding to the top of the core.
    i_depth : integer
        The index of the radius that is closest to radius[n] - depth. This is
        returned only when the optional parameter depth is specified.

    Parameters
    ----------
    filename : str
        Name of the input file
    depth : float, optional, default = None
        If specified, the returned value i_depth will correspond to the index
        of the radius vector that is closest to radius[n] - depth.
    quiet : bool, optional, default = True
        If False, print summary information about the model.

    Notes
    -----
    This routine reads in a "deck" reference interior model for use in the
    routine HydrostaticShapeLith. Note that the latter routine requires the
    density to be a specified constant value between two radii (see Figure 1 of
    Wieczorek et al. 2019), whereas the format of the "deck" files provides the
    density at each specified radius. Thus, this routine sets density[i] equal
    to the average of the density at radius[i] and radius[i+1] in the original
    file.

    The deck file must be formatted as follows. The first line is a header
    line, which is output only when quiet is set to False. The second line
    lists three parameters: (1) ifanis, which is 0 or 1 for an isotropic or
    anisotropic model, respectively, (2) tref, which is the reference period in
    seconds of the model for the physical dispersion correction, and (3)
    ifdeck, which is 1 for tabular data or 0 for a polynomial model
    (unsupported). The third line lists the following: (1) N, the number of
    model entries, excluding the first three lines of the file, (2) the line
    number, starting from 1, corresponding to the top of the solid inner core,
    (3) the line number corresponding to the top of the fluid core, and (4) the
    line number corresponding to the top of the mantle. Following this, each
    line must contain the radius and density (other values are not used).
    Discontinuities at the inner-outer core boundary, core-mantle interface and
    crust-mantle interface are represented by two adjacent lines with identical
    radii providing the density on either side of the discontinuity.
    '''

    with open(filename, 'r') as f:
        lines = f.readlines()
        data = lines[1].split()
        if float(data[2]) != 1:
            raise RuntimeError('Program not capable of reading polynomial ' +
                               'files')
        num_file = int(lines[2].split()[0])
        crust_index_file = int(lines[2].split()[3])
        core_index_file = int(lines[2].split()[2])
        i_crust_file = crust_index_file - 1
        i_core_file = core_index_file - 1

        radius = np.zeros(num_file)
        rho = np.zeros(num_file)
        num = 0

        for i in range(0, num_file-1):
            data = lines[i+3].split()
            rb = float(data[0])
            rhob = float(data[1])
            data = lines[i+4].split()
            rt = float(data[0])
            rhot = float(data[1])

            if rb == rt:
                if i == i_core_file:
                    i_core = num
                if i == i_crust_file:
                    i_crust = num
            else:
                radius[num] = rb
                rho[num] = (rhot + rhob) / 2.
                num += 1

        radius[num] = rt
        rho[num] = 0.  # the density above the surface is zero
        num += 1
        n = num - 1
        radius = radius[:n+1]
        rho = rho[:n+1]

        rho_mantle = rho[i_crust-1]
        rho_core = rho[i_core-1]
        r0_model = radius[n]

        if depth is not None:

            for i in range(0, n+1):
                if radius[i] <= (r0_model - depth) and \
                        radius[i+1] > (r0_model - depth):
                    if radius[i] == (r0_model - depth):
                        i_depth = i
                    elif (r0_model - depth) - radius[i] <= radius[i+1] -\
                            (r0_model - depth):
                        i_depth = i
                    else:
                        i_depth = i + 1
                    break

        if quiet is False:
            print(lines[0].strip())
            print('Surface radius of model (km) = {:f}'
                  .format(r0_model / 1.e3))
            print('Mantle density (kg/m3) = {:f}'.format(rho_mantle))
            print('Mantle radius (km) = {:f}'.format(radius[i_crust]/1.e3))
            print('Core density (kg/m3) = {:f}'.format(rho_core))
            print('Core radius (km) = {:f}'.format(radius[i_core]/1.e3))

            if depth is not None:
                print('Assumed depth of optional interface (km) = {:f}'
                      .format(depth / 1.e3))
                print('Actual depth of optional interface in discretized '
                      'model (km) = {:f}'.format((r0_model - radius[i_depth])
                                                 / 1.e3))

        if depth is not None:
            return radius, rho, i_crust, i_core, i_depth
        else:
            return radius, rho, i_crust, i_core
