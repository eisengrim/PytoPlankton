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
import matplotlib.tri as Tri
from astropy import Time
import pyproj

# local imports
from temporal_utilities import *
from spatial_utilities import *

class _load_settings:
    """
    Loads the settings for the Plankton object, and stores it in a subclass.
    It is necessary to change these parameters before running the Lagrangian
    tracker code.
    """
    def __init__(self, grid_file, locs_file, kwargs):
        params_default = {lps : 60,
                          ops : 60,
                          start : 10,
                          total : 30,
                          solver : 'rk4',
                          interp : 'linear',
                          lon : [],
                          lat : []}

        # update settings attributes
        for (param, default) in params_default.iteritems():
            setattr(self, param, kwargs.get(param, default))

        self.locs_file = locs_file
        self.grid_file = grid_file
                    

class _load_fvcom_grid:
    """
    Loads the grid data for the subclass Grid. Calculates information about
    the grid needed to increase particle tracking speed. Data is from FVCOM.

    --VARIABLES--
    x : x-coordinte at nodes (m); (nnode)
    y : y-coordinate at nodes (m); (nnode)
    xc : x-coordinate at element (m); (nele)
    yc : y-coordinate at element (m); (nele)
    lon : longitude at node (deg); (nnode)
    lat : latitude at node (deg); (nnode)
    lonc : longitude at element (deg); (nele)
    latc : latitude at element (deg); (nele)
    h : bathymetric height (m); (ntime, nnode)
    hele : bathymetric heights at elements (m); (ntime, nnode)
    siglay : sigma layers; (nlevel, nnode)
    siglev : sigma levels; (nlevel+1, nnode)
    nbn : nearest bounding nodes, surrounding indices; (3, nnodes) [nv]
    nbe : nearest bounding elements, surrounding indices; (3, nele)
    nnode : number of nodes, integer
    nele : number of elements, integer
    nlevel : number of vertical elements, integer
    ntime: number of time elements, integer
    nsiglay : siglay dimension, integer
    nsiglev : siglev dimension, integer
    zlay : depth of each node at each sigma layer (m)
    zlev : depth of each node at each sigma level (m)
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
    def __init__(self, gridpath, settings, debug=False):
        # look for directory and ncfile
        if debug:
            print 'retrieving data from ' + gridpath + '...'
        if not osp.exists(gridpath):
            sys.exit('...the file {} was not found.'.format(gridpath))

        try:
            data = sio.netcdf.netcdf_file(gridpath, 'r', mmap=True)
        except:
            data = nc.Dataset(gridpath, 'r', format='NETCDF4_CLASSIC')

        if debug:
            print '\tloading grid variables...'

        # load grid variables from raw data
        # add __slots__?
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
        if debug:
            print '\tloading layers...'
            
        lsig, nodes = self.siglev.shape
        self.zlay = siglay[0:lsig-1, 0:1]
        self.zlev = siglev[0:lsig, 0:1]
        self.nlevel = self.siglay.shape[0]
        self.nele = self.lonc.shape[0]
        self.nnode = self.lon.shape[0]
        self.nsiglay = len(self.zlay)
        self.nsiglev = len(self.zlev)

        # add to code, add to documentation
        # load rest of grid variables here
        # self.n_isonb = np.zeros(self.nnode)
        # self.e_isonb = np.zeros(self.nele)
        # self.uin = np.zeros(self.nele+1, self.nsiglay)
        # self.vin = np.zeros(self.nele+1, self.nsiglay)
        # self.win = np.zeros(self.nele+1, self.nsiglay)
        # self.hin = np.zeros(self.nnode)
        
        # self.unc1, wnc2, uin, vin, hin, win, etc??

        # set flag for boundary nodes and elements
        # self.e_isonb[np.where(self.nbe == 0)] == 1
        # self.n_isonb[np.where(self.nbn[np.where(self.nbe

        # calculate global minimums and maximums
        if debug:
            print '\tcalculating global minimums and maximums...'
            
        self.xmin = np.min(self.x)
        self.xmax = np.max(self.x)
        self.ymax = np.max(self.y)
        self.ymin = np.min(self.y)

        # shift grid to upper right Cartesian
        self.x = self.x - self.xmin
        self.y = self.y - self.ymin

        # load time information, julian and matlab times
        try:
            self.mjd = data.variables['time'].data
        except AttributeError:
            # exception for nc.Dataset
            self.mjd = data.variables['time']
        
        if debug:
            print '\t...passed'


class _load_mat_grid:
    """
    Loads the grid data from MATLAB scatter data.

    --VARIABLES--
    lon : 
    lat :
    mlon :
    mlat :
    time :
    uttc :
    vttc :
    uuss :
    vuss : 
    uwnd :
    vwnd :
    proj :
    x :
    y :
    mx :
    my :
    mask :
    ustep1 :
    vstep1 :
    ustep2 :
    vstep2 :
    
    """
    def __init__(self, gridpath, settings, debug=False):
        # look for directory and matlab file
        if debug:
            print 'retrieving data from ' + gridpath + '...'
        if not osp.exists(gridpath):
            sys.exit('...the file {} was not found.'.format(gridpath))
        # enclose in try statement?
        data = sio.loadmat(gridpath)

        if debug:
            print '\tloading grid variables...'
            
        # load grid variables
        # add __slots__?
        gridvars = ['lon', 'lat', 'time', 'uttc', 'uuss', 'uwnd', 'vttc', \
                    'vuss', 'vwnd']

        for key in gridvars:
            setattr(self, data[key])
        
        # define the lcc projection
        if debug:
            print '\tsetting lcc projection...'
            
        self._xmax = np.nanmax(self.lon)
        self._xmin = np.nanmin(self.lon)
        self._ymax = np.nanmax(self.lat)
        self._ymin = np.nanmax(self.lat)

        self._xavg = (self._xmax + self._xmin) * 0.5
        self._yavg = (self._ymax + self._ymin) * 0.5
        self._ylower = (self._ymax - self._ymin) * 0.25 + self.ymin
        self._yupper = (self._ymax - self._ymin) * 0.75 + self.ymin

        self._projstr = 'lcc +lon_0=' + str(self._xavg) + ' +lat_0=' \
                       + str(self._yavg) + ' +lat_1=' + str(self._ylower) \
                       + ' +lat_2=' + str(self._yupper)
        self.proj = pyproj.Proj(proj=self._projstr)

        self.mlon, self.mlat = np.meshgrid(self.lon, self.lat)
        self.lon, self.lat = self.mlon.flatten(), self.mlat.flatten()
        self.x, self.y = self.proj(self.lon, self.lat)
        self.mx, self.my = self.proj(self.mlon, self.mlat)

        if debug:        
            print '\tflattening data fields...'

        # flatten the data as it is basically scattered in xy
        # self.uwnd = flattime(self.uwnd)
        # self.vwnd = flattime(self.vwnd)
        # self.uuss = flattime(self.uuss)
        # self.vuss = flattime(self.vuss)
            
        self.uttc = flattime(self.uttc)
        self.vttc = flattime(self.vttc)

        # save the mask
        self.mask = np.isnan(self.uttc)[:,:,0]

        self.ustep1 = self.uttcf[0,:]
        self.vstep1 = self.vttcf[0,:]
        self.ustep2 = self.uttcf[0,:]
        self.vstep2 = self.vttcf[0,:]

        self.mjd = dn2mjd(self.time)
        
        if debug:
            print '\t...passed.'

        
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
    def __init__(self, grid, settings, debug=False):
        if debug:
            print '\tloading time variables...'

        # Ensures internal step and output step are scalar multiples
        self.internalstep = np.float64(settings.lps)
        self.outputstep = np.float64(settings.ops)
        if not np.mod(internalstep, outputstep) == 0:
            sys.exit('output and internal steps must be evenly divisible.')

        self.mjd = grid.mjd
        
        self.totalsteps = np.float64(settings.total)
        self.startstep = np.float64(settings.start)

        if debug:
            print '\tconverting time measurements...'
        self.mdn = mjd2dn(self.mjd)
        self.ntime = self.mjd.shape[0]

        # convert mjd to jd, then to a datetime
        self.jd = mjd2jd(self.mjd)
        self._dates = jd2dt(self.jd)

        if debug:
            print '\tcalculating stepping...'
            
        # calculate stepping of input data, internal and output stepping
        self.instp = np.round((self.jd[1]-self.jd[0]).sec)
        self.dti = self.instp/self.internalstep
        self.dtout = self.instp/self.outputstep

        self.finishstep = self.startstep + self.totalsteps

        # add self.itout, iint, i2, int2, outt?

        # number of outputs to complete
        self.nouts = (self.totalsteps * self.outputstep) + 1

        if debug:
            print '\t...passed'


class _load_part:
    """
    Loads the Lagrangian Particle subclass.

    --VARIABLES--
    """
    def __init__(self, grid, time, settings, debug=False):
        if debug:
            print '\tloading particle variables...'

        locs_file = settings.locs_file

        if debug:
            print '\tretrieving initial positions from ' + locs_file + '...'
        if not osp.exists(locs_file):
            sys.exit('cannot find location file {}.'.format(locs_file))

        # set up initial positions, ignore missing data
        locs = np.genfromtxt(locs_file, comments='#', autostrip=True)
        
        if not len(settings.lon) == len(settings.lat):
            sys.exit('number of longitudes and latitudes given must match.')

        # load lon/lat and elem given from settings
        if settings.lon and settings.lat:
            lon = np.asarray(settings.lon)
            lat = np.asarray(settings.lat)
            locs_opt = np.vstack((lon,lat)).T
            self.locs = np.vstack((locs,locs_opt))        

        # initialize arrays for lagrangian particle(s)
        self.locs = self.locs[~np.isnan(self.locs).any(axis=1)]
        self.nparts = self.locs.shape[0]

        n = self.nparts
        self.x = np.zeros(n, time.nouts)
        self.y = np.zeros(n, time.nouts)

        
    # initialize arrays for lagrangian particle, set particle positions
    # lag.nparts = n
    # lag.time = np.zeros(time.nouts)
    # lag.x = np.zeros(n, # number times out?
    # lag.y
    # lag.z
    # lag.u
    # lag.w
    # lag.v
    # lag.sig_pos
    # lag.host
    # lag.indomain
    # lag.h
    # lag.free_sfc

    # difference between position / velocity and their absolutes?

