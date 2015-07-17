#!/usr/bin/env python
# encoding: utf-8

# library imports
from __future__ import division
import sys
import os
import os.path as osp
from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import scipy.io as sio
import scipy.interpolate as si
import matplotlib.tri as Tri
import seaborn as sns


# local imports
from temporal_utilities import *


def find_index(lon, lat, grid):
    """
    Find the index of a long / lat combination.
    """
    # take it from PySeidon


def interp(settings, grid, infield, lon, lat, trinodes, debug=False):
    """
    Interpolation of any variable at any position...

    --VARIABLES--
    settings : Plankton Settings subclass
    grid : Plankton Grid subclass
    infield : grid variable to be interpolated
    lon : longitude to interpolate
    lat : latitiude to interpolate
    trinodes : FVCOM triangle nodes, numpy array (3, nele)
    """

    if debug:
        print '\tinterpolating at point...'
        
    # take from PySeidon
    if settings.interp == 'gridded':
        outvar = _gridded(infield, grid, lon, lat, debug=False)

    if settings.interp == 'linear':
        outvar = _linear(infield, grid, lon, lat, trinodes, debug=False)


def _gridded(var, grid, lon, lat):
    """
    Gridded interpolation.
    """
    return si.griddata(np.vstack([grid.lon, grid.lat]).T, \
                       var, np.vstack([lon, lat]).T)


def _linear(var, grid, lon, lat, trinodes):
    """
    Linear interpolation of variable at any location from PySeidon
    """
    
    triInd = trinodes
