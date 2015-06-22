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

    --VARIABLES--
    x : x-coordinte at nodes (m); (ntime, nnode)
    y : y-coordinate at nodes (m); (ntime, nnode)
    xc : x-coordinate at element (m); (ntime, nele)
    yc : y-coordinate at element (m); (ntime, nele)
    lon : longitude at node (deg); (ntime, nnode)
    lat : latitude at node (deg); (ntime, nnode)
    lonc : longitude at element (deg); (ntime, nele)
    latc : latitude at element (deg); (ntime, nele)
    h : bathymetric height (m); (ntime, nnode)
    hele : bathymetric heights at elements (m); (ntime, nnode)
    siglay : sigma layers; (nlevel, nnode)
    siglev : sigma levels; (nlevel+1, nnode)
    nbn : nearest bounding nodes, surrounding indices; (3, nele) [nv]
    nbe : nearest bounding elements, surrounding indices; (3, nele)
    nnode : number of nodes, integer
    nele : number of elements, integer
    nlevel : number of vertical elements, integer
    ntime: number of time elements, integer
    nsiglay : siglay dimension, integer
    nsiglev : siglev dimension, integer
    zlay :
    zlev :
    a1u, a2u, awx, awy, aw0 : grid parameters

    flag isonb set if node / element is on a boundary:
        n_isonb(i) = 0 : node in the interior computational domain
        n_isonb(i) = 1 : node is on the solid boundary
        n_isonb(i) = 2 : node is on the open boundary
        e_isonb(i) = 0 : element is in the interior computational domain
        e_isonb(i) = 1 : element is on the solid boundary
        e_isonb(i) = 2 : element is on the open boundary
        e_isonb(i) = 3 : element with 2 solid boundary edges
    """
    def __init__(self, data, settings, debug=False):
        if debug:
            print '\tloading grid variables...'

        # load grid variables from raw data
        gridvars = ['x', 'y', 'xc', 'yc', 'h', 'a1u', 'a2u', 'awx', 'awy', \
                    'aw0', 'siglev', 'siglay']
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
        self.n_isonb = np.zeros(self.nnode)
        self.e_isonb = np.zeros(self.nele)
        self.uin = np.zeros(self.nele+1, self.nsiglay)
        self.vin = np.zeros(self.nele+1, self.nsiglay)
        self.win = np.zeros(self.nele+1, self.nsiglay)
        self.hin = np.zeros(self.nnode)

        # self.unc1, wnc2? 

        # calculate global minimums and maximums
        self.xmin = np.amin(self.x)
        self.xmax = np.amax(self.x)
        self.ymax = np.amax(self.y)
        self.ymin = np.amin(self.y)

        # shift grid to upper right Cartesian
        self.x = self.x - self.xmin
        self.y = self.y - self.ymin

        # calculate approximate element coordinates
        self.xc = (self.x[self.nbn[:,1]] + self.x[self.nbn[:,2]] \
                   + self.x[self.nbn[:,3]]) / 3
        self.yc = (self.y[self.nbn[:,1]] + self.y[self.nbn[:,2]] \
                   + self.y[self.nbn[:,3]]) / 3

        # set flag for boundary nodes and elements
        # self.e_isonb[np.where(self.nbe == 0)] == 1
        
        # find the closest ## elements to each element and write-save them?
        
        if debug:
            print '\t...passed'

        
class _load_time_var:
    """
    Loads the time variables subclass Time.

    --VARIABLES--
    internalstep : number of loops to go through each fvcom time step
    outputstep : number of times to output each fvcom time step
    totalsteps : total tracking time in fvcom time steps
    startstep : starting step of the lagrangian tracker
    mjd : time in modified julian date
    jd : time in julian date
    mdn : time expressed as a matlab datenum
    instp : stepping of input data from fvcom
    dti : internal stepping 
    dtout : output stepping
    ntime : dimension of time variable
    finishstep : last step of lagrangian model run
    nouts : total number of outputs
    """
    def __init__(self, data, settings, debug=False):
        if debug:
            print '\tloading time variables...'

        # Ensures internal step and output step are scalar multiples
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

        # number of outputs to complete
        self.nouts = (self.totalsteps * self.outputstep) + 1
        
        if debug:
            print '\t...passed'
                                       
        
class _load_part:
    """
    Loads the Lagrangian Particle subclass.

    --VARIABLES--
    """
    def __init__(self, data, time, settings, debug=False):
    
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
    lag.npts = settings['number_of_particles']
    lag.time = np.zeros(time.nouts)
    
    
