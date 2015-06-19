#!/usr/bin/python2.7
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
import matplotlib.tri as Tri
from astropy import Time

# local imports
from temporal_utilities import *
from spatial_utilities import *

class _load_settings:
    """
    Loads the settings for the Plankton object, and stores it in a subclass.
    It is necessary to change these parameters before running the Lagrangian
    tracker code. This is done by simply calling the Settings subclass and
    assigning a new value to the key.
    """
    def __init__(self, **kwargs):
        # set defaults, update with given kwargs
        params_default = {'loops_per_step' : 60,
                          'output_per_step' : 60,
                          'start_step' : 4320,
                          'number_of_steps' : 30,
                          'number_of_particles' : 10,
                          'sigstart_option' : 4,
                          'x_box_size' : 100,
                          'y_box_size' : 100,
                          'radius' : 10,
                          'pos_option': 1,
                          'element_loc' : 41577,
                          'x_loc' : None,
                          'y_loc' : None, 
                          'lon' : None,
                          'lat' : None,
                          'bound_box' = False,
                          'bounds' = []}

        for (param, default) in params_default.iteritems():
            setattr(self, param, kwargs.get(param, default))

    
class _load_grid:
    """
    Loads the grid data for the subclass Grid. Calculates information about
    the grid needed to increase particle tracking speed.

    Description of variables:
    
    """
    def __init__(self, data, settings, debug=False):
        if debug:
            print '\tloading grid variables...'

        # load grid variables from raw data
        gridvars = ['x', 'y', 'h', 'a1u', 'a2u', 'awx', 'awy', 'aw0', \
                    'siglev', 'siglay']
        for key in gridvars:
            try:
                setattr(self, key, data.variables[key].data)
            except AttributeError:
                # exception for nc.Dataset type
                setattr(self, key, data.variables[key])

        # special treatment for triele and trinodes
        datavar = data.variables.keys()
        if 'trinodes' in datavar:
            try:
                setattr(self, 'nbn', data.variables['trinodes'].data)
            except AttributeError:
                setattr(self, 'nbn', data.variables['trinodes'])
        else:
            try:
                self.nbn = np.transpose(data.variables['nv'].data) - 1
            except AttributeError:
                self.nbn = np.transpose(data.variables['nv'].data) - 1
        if 'triele' in datavar:
            try:
                setattr(self, 'nbe', data.variables['triele'].data)
            except AttributeError:
                setattr(self, 'nbe', data.variables['triele'])
        else:
            try:
                self.nbe = np.transpose(data.variables['nbe'].data) - 1
            except AttributeError:
                self.nbe = np.transpose(data.variables['nbe'].data) - 1

        # approximate bathymetry at the elements
        self.hele = (self.h[self.nbn[0,:]-1] + self.h[self.nbn[1,:]-1] \
                     + self.h[self.nbn[2,:]-1]) / 3
        
        # load information relating to sigma layers and # of nodes/elements
        lsig, nodes = self.siglev.shape
        self.zlay = siglay[0:lsig-1, 0:1]
        self.zlev = siglev[0:lsig, 0:1]
        self.nlevel = self.siglay.shape[0]
        self.nele = self.lonc.shape[0]
        self.nnode = self.lon.shape[0]
        self.nsiglay = len(self.zlay)
        self.nsiglev = len(self.zlev)

        # load rest of grid variables here
        self.isonb = np.zeros(self.nele)
        self.isbce = np.zeros(self.nele)
        self.uin = np.zeros(self.nele+1, self.nsiglay)
        self.vin = np.zeros(self.nele+1, self.nsiglay)
        self.win = np.zeros(self.nele+1, self.nsiglay)
        self.hin = np.zeros(self.nnode)

        # self.unc1, wnc2?
        # triangleedgegrid
        
        if debug:
            print '\t...passed'

        
class _load_time_var:
    """
    Loads the time variables subclass Time.

    Description of variables:
    
    """
    def __init__(self, data, settings, debug=False):
        if debug:
            print '\tloading time variables...'

        # ensures internal step and output step are scalar multiples
        self.internalstep = np.float64(settings.loops_per_step)
        self.outputstep = np.float64(settings.output_per_step)
        if not np.mod(instp, otstp) == 0:
            sys.exit('output and internal steps must be evenly divisible.')

        self.totalsteps = np.float64(settings.number_of_steps)
        self.startstep = np.float64(settings.start_step)

        # load time information, julian and matlab times
        try:
            self.mjd = data.variables['time'].data
        except AttributeError:
            # exception for nc.Dataset
            self.mjd = data.variables['time']

        self.mdn = mjd2dn(self.mjd)
        self.ntime = self.mjd.shape[0]

        # convert mjd to jd, then to a datetime
        self.jd = mjd2jd(self.mjd)
        self._dates = jd2dt(self.jd)

        # calculate stepping of input data, internal and output stepping
        self.instp = np.round((self.jd[1]-self.jd[0]).sec)
        self.dti = self.instp/self.internalstep
        self.dtout = self.instp/self.outputstep

        self.finishstep = self.startstep + self.totalsteps

        # self.itout, iint, i2, int2, outt?
        
        if debug:
            print '\t...passed'
        
        
class _load_part:
    """
    Loads the Lagrangian Particle subclass.

    """
    def __init__(self, data, settings, debug=False):
    
    if debug:
        print '\tloading particle variables...'

    # place particles in initial positions based on selections made
    if settings.sigstart_option == 1:
        pass
    elif settings.sigstart_option == 2:
        pass
    elif settings.sigstart_option == 3:
        pass
    elif settings.sigstart_option == 4:
        pass

    # initialize arrays for lagrangian particle, set particle positions


    
