#!/usr/bin/env python
# encoding: utf-8

# library imports
from __future__ import division
import sys
import os
import os.path as osp
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4 as nc
import scipy as sp
import scipy.io as sio
import scipy.interpolate as itp
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as tic
import seaborn as sns
# import argparse as arp
# import cPickle as cpk
# import pickle as pkl
# import gc
import pyproj

# local imports
from planktonIO import *

class Plankton:
    """
    Creates a Plankton tracker object in PytoPlankton for particle tracking.

    ---STRUCTURE---
              /-Grid : Holds all grid data.
             |--Settings : Holds all specified settings.
    Plankton_|--Time : Holds time-related data and settings for particles.
             |--Particle : Holds various data for particles.
              \-Write : Save options.

    ---INPUT---
    - gridpath = Path to the grid file that holds the grid data. Either a
        netCDF4 file from FVCOM, a pickle file or a MATLAB file.
    - locspath = Path to the .dat file that holds all of the initial
        positions of the particles in lonlat or elem format. If in lonlat,
        two columns of lon and lat are separated by a single space each and 
        each particle has its own row. If in element format, elements are
        given in a single column with each row dedicated to a particle.
        Note that this may be augmented by passing the parameters lonlat 
        or elem.
        
    ** Optional Parameters:
    - lps : Loops per step. Number of loops to be completed in one FVCOM time
        step, or the linear time interpolation ratio between nc steps.
        Default value is 60, so if the FVCOM model saves data every 1 hour,
        particle positions are calculated every minute. Integer.
    - ops : Output per step. Number of times to save data in one FVCOM time
        step, similar in principle to lps. Note that in order to save, it 
        is required that mod(lps, ops)==0. Default value is 60. Integer.
    - start : The starting step of the tracker. Default is 30 time steps. 
    - total : The total number of steps for which to track the particles.
        The default value is 30 time steps. Integer.
    - lon / lat : initial positions given as lon / lat lists.
    - solver : Current options for the numerical solver are: 
             'fwd_euler' (default): for Forward Euler method.
             'rk4' : perform a Runge-Kutta 4 Method.
    - interp_method : Interpolation method. Options include:
             'linear' (default) : perform a linear interpolation.
             'gridded' : perform an interpolation of uniformly gridded data.
    - debug = Boolean. Increase output verbosity if True.
    
    * Coordinates are given as decimal degrees East and North.
    * Times are printed out in compliance with ISO 8601 format standards,
      <date>T<time>, with T as the delimiter: YYYY-MM-DDTHH:MM:SS.ss (UTC).

    """

    def __init__(self, gridpath, locspath, debug=False, **kwargs):
        '''
        Initializes Plankton class.
        '''
        self._debug = debug
        if debug:
            print '-debug mode turned on-'

        # declares settings to be used...
        if debug:
            print 'loading settings...'
        params = kwargs
        self.Settings = _load_settings(gridpath, locspath, params)

        # load file (add opendap functionality?)
        if debug:
            print 'looking for grid data...'

        # load grid. note this may be a large file, no need to move to input
        if gridpath.endswith('.nc'):
            self.Grid = _load_fvcom_grid(gridpath,
                                         self.Settings,
                                         debug=self._debug)

        elif gridpath.endswith('.mat'):
            self.Grid = _load_mat_grid(gridpath,
                                       self.Settings,
                                       debug=self._debug)

        else:
            sys.exit('...file format for grid data is not supported.')
            
        # load data into structures
        if debug:
            print 'initializing plankton object...'
        try:           
            # load time data and initialize particle
            self.Time = _load_time_var(self.Grid,
                                       self.Settings,
                                       debug=self._debug)
            self.Particle = _load_part(self.Grid,
                                        self.Time,
                                        self.Settings,
                                        debug=self._debug)

        except MemoryError:
            sys.exit('...data too large for machine memory.')

        # add plot and functions here from local import?
        # special methods like __del__, __new__, __add__ here?
