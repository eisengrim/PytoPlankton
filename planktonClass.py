#!/usr/bin/python2.7
# encoding: utf-8

# library imports
from __future__ import division
import sys
import os
import argparse as arp
import os.path as osp
from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import scipy as sp
import scipy.io as sio
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import matplotlib.ticker as tic
import seaborn as sns
import cPickle as cpk
import pickle as pkl
import gc

# local imports
from setupPlankton import *


class Plankton:
    """
    Creates a Plankton tracker object in PytoPlankton for particle tracking.

    ---STRUCTURE---
              /-Grid : Holds all information and data related to the grid.
             |--Settings : Holds various settings for initialization.
    Plankton_|--Time : Holds time-related data and settings for particles.
             |--Particle : Holds various spatial data for particles.
             |--Plot : Useful functions for plotting.
              \-Write : Functions to save the data.

    ---INPUT---
     - filename : Path to input file to be read in. Contains
        all relevant settings for initializing the object in netCDF4 format.
     - kwargs : Any optional parameter can be changed here or later by
        calling Plankton.Settings.
     - debug = [ True | False ] : Increases output verbosity if True.

    ---PARAMETERS---
    Additional keyword arguments are defined as follows:
     - loops_per_step : Number of loops to be completed in one FVCOM time
         step, or the linear time interpolation ratio between nc steps.
         Default value is 60, so if the FVCOM model saves data every 1 hour,
         particle positions are calculated every minute. If FVCOM data is
         produced every minute, than the calculations are done for every
         second.
     - output_per_step : Number of times to save the data in one FVCOM time
         step, similar in principle to loops_per_step. Note that in order to
         save, it is required that mod(loops_per_step, output_per_step)==0.
         Default value is 30, or every other second if the FVCOM output is in
         minutes.
     - start_step : The starting step of the tracker. Default is 4320 time
         steps, equal to 3 days in minutes.
     - number_of_steps : The total number of steps to track the particle for.
         The default value is 30, which is 30 minutes if the FVCOM output is
         every minute.
     - number_of_particles : Number of particles to track.
     - sigstart_option : The option to use for placing the particles in
         their initial positions. Defualt is 4.
             1 : A group of particles is randomly placed in a specified box
                 around a location, given as a position.
             2 : A group of particles is randomly placed within a circle of
                 a given radius around a location, given as a position.
             3 : Start particles at the centres of the elements found within
                 a specified box of a certain size around a position.
             4 : The particles begin at a specific position given.
     - x_box_size, y_box_size : The size of the box in which the particles
         are generated. Units are in metres. Default is 100 m.
     - radius : The radius of the circle in which the particles are generated
         in metres. Default is 10 m.
     - pos_option : The option to use for defining locations. Default is 1.
             1 : Positions are given as element numbers.
             2 : Positions are given in terms of Cartesian coordinates.
             3 : Positions are given in terms of geographic coordinates*.
     - element_loc : The element used in sigstart_option, either as the
         centre of a box, circle, or an independent point. Default is 41577.
         This option is only used if positions == 1.
     - x_loc, y_loc : The x and y coordinates to be used in the options above,
         if positions = 2. Defaults are None and None.
     - lon, lat : The longitudal and latitudal coordinates to be used in the
         options above, if positions == 3. Default is None and None.
     - bound_box : Boolean. If true, the box will be bounded not by size, but
         in terms of longitudes and latitudes. Default is False.
     - bounds : The longitudes and latitudes to be used in creating the box.
         The box is given in the form [min_long, max_long, min_lat, max_lat].
         This list is by default empty.

    * Coordinates are given as decimal degrees East and North.
    * Times are printed out in compliance with ISO 8601 format standards,
      <date>T<time>, with T as the delimiter: YYYY-MM-DDTHH:MM:SS.ss (UTC).

    """

    def __init__(self, filename, debug=False, **kwargs):
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
        self.Settings = _load_settings(params)
        
        # load netcdf file (add pickle and opendap functionality?)
        if not filename.endswith('.nc'):
            sys.exit('...file format is not supported.')

        # look for directory and ncfile
        if debug:
            print 'retrieving data from ' + filename + '...'
        if not osp.exists(filename) or not osp.isfile(filename):
            sys.exit('...the file {} was not found.'.format(filename))
        try:
            self._data = sio.netcdf.netcdf_file(filename, 'r', mmap=True)
        except:
            self._data = nc.Dataset(filename, 'r', format='NETCDF4_CLASSIC')

        # load data into structures
        if debug:
            print 'initializing plankton object...'
        try:
            self.Grid = _load_grid(self._data,
                                   self.Settings,
                                   debug=self._debug)
            self.Time = _load_time_var(self._data,
                                       self.Settings,
                                       debug=self._debug)
            self.Particle = _load_part(self._data,
                                        self.Time,
                                        self.Settings,
                                        debug=self._debug)

        except MemoryError:
            sys.exit('...data too large for machine memory.')

        # add plot and functions here from local import?
        # special methods like __del__, __new__, __add__ here?
