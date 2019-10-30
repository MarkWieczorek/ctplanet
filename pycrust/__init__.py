"""
pycrust
=======

pycrust provides several functions for generating crustal thickness maps of a
planet from gravity and topography data, and the calculation of hydrostatic
relief along density interfaces beneath the lithosphere.

Contents

    pyMoho                    Calculate relief using a constant density crust
                              and mantle.
    pyMohoRho                 Calculate relief using a constant density mantle
                              and a variable density crust.
    HydrostaticShapeLith      Calculate the relief of hydrostatic interfaces
                              beneath the lithosphere along with the predicted
                              gravity, taking into account rotation and/or
                              tides using the approach of Wieczorek et al.
                              (2019).
    HydrostaticShape          Calculate the relief of hydrostatic interfaces
                              and predicted gravity of a rotating hydrostatic
                              planet using the approach of Wieczorek et al.
                              (2019).
    InertiaTensor_from_shape  Calculate the inertia tensor given a radial
                              density profile and shape of each interface.
    InertiaTensor_from_C      Calculate the inertia tensor given the polar
                              moment of inertia and the gravitational potential
                              coefficients.
    moi                       Calculate the mean, normalized, moment of inertia
                              up to index n.
    ReadRefModel              Read the reference interior model file.

"""
__version__ = '0.1'
__author__ = 'Mark Wieczorek'


from .pyMoho import pyMoho
from .pyMoho import pyMohoRho

from .Hydrostatic import HydrostaticShapeLith
from .Hydrostatic import HydrostaticShape

from .InertiaTensor import InertiaTensor_from_shape
from .InertiaTensor import InertiaTensor_from_C
from .InertiaTensor import moi

from .ReadRefModel import ReadRefModel


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['pyMoho', 'pyMohoRho', 'HydrostaticShapeLith', 'HydrostaticShape',
           'InertiaTensor_from_shape', 'InertiaTensor_from_C', 'moi',
           'ReadRefModel']
