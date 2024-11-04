"""
ctplanet
=======

ctplanet provides several functions for generating crustal thickness maps of a
planet from gravity and topography data, and the calculation of hydrostatic
relief along density interfaces beneath the lithosphere.

Notes

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
from importlib.metadata import version, PackageNotFoundError

from .Moho import pyMoho
from .Moho import pyMohoRho

from .Hydrostatic import HydrostaticShapeLith
from .Hydrostatic import HydrostaticShape

from .InertiaTensor import InertiaTensor_from_shape
from .InertiaTensor import InertiaTensor_from_C
from .InertiaTensor import moi

from .ReadRefModel import ReadRefModel

del Moho  # noqa: F821
del Hydrostatic  # noqa: F821
del InertiaTensor  # noqa: F821

try:
    __version__ = version('ctplanet')
except PackageNotFoundError:
    # package is not installed
    pass

__author__ = 'Mark Wieczorek'

__all__ = ['pyMoho', 'pyMohoRho', 'HydrostaticShapeLith', 'HydrostaticShape',
           'InertiaTensor_from_shape', 'InertiaTensor_from_C', 'moi',
           'ReadRefModel']
