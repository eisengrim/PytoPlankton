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
import time
import pyproj
from astropy.time import Time

# local imports
from temporal_utilities import *
from spatial_utilities import *

class load_settings:
    """
    Loads the settings for the Plankton object, and stores it in a subclass.
    It is necessary to change these parameters before running the Lagrangian
    tracker code.
    """
    def __init__(self, grid_file, locs_file, out_file, kwargs):
        params_default = {'lps' : 60,
                          'ops' : 60,
                          'start' : 10,
                          'total' : 30,
                          'solver' : 'rk4',
                          'interp' : 'linear',
                          'lon' : [],
                          'lat' : []}

        # update settings attributes
        for (param, default) in params_default.iteritems():
            setattr(self, param, kwargs.get(param, default))

        self.locs_file = locs_file
        self.grid_file = grid_file
        self.out_path = out_file
        self.finish = self.start + self.total


class load_fvcom_grid:
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
    trinodes : nearest bounding tri nodes, surrounding indices; (3, nele) [nv]
    triele : nearest tri elements, surrounding indices; (3, nele) [nbe]
    nnode : number of nodes, integer
    nele : number of elements, integer
    nlevel : number of vertical elements, integer
    ntime: number of time elements, integer
    nsiglay : siglay dimension, integer
    nsiglev : siglev dimension, integer
    zlay : depth of each node at each sigma layer (m)
    zlev : depth of each node at each sigma level (m)

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
            print 'retrieving data from ' + osp.basename(gridpath) + '...'
        if not osp.exists(gridpath):
            print '...the file {} was not found.'
            sys.exit()

        try:
            data = sio.netcdf.netcdf_file(gridpath, 'r', mmap=True)
        except:
            data = nc.Dataset(gridpath, 'r', format='NETCDF4_CLASSIC')

        if debug:
            print 'loading grid variables...'

        # load grid variables from raw data
        # add __slots__?
        datavar = data.variables.keys()
        # determine if 3d
        if 'h' in datavar:
            self.threeD = True
        
        gridvars = ['x', 'y', 'xc', 'yc', 'lon', 'lat', 'lonc', 'latc'] 
        if self.threeD:
            gridvars = gridvars + ['h', 'siglev', 'siglay']

        if debug:
            print '\tsettings attributes...'
        for key in gridvars:
            try:
                setattr(self, key, data.variables[key].data)
            except AttributeError:
                # exception for nc.Dataset type
                setattr(self, key, data.variables[key])

        # special treatment for triele and trinodes
        if 'trinodes' in datavar:
            try:
                setattr(self, 'trinodes', data.variables['trinodes'].data)
            except AttributeError:
                setattr(self, 'trinodes', data.variables['trinodes'])
        else:
            try:
                self.trinodes = np.transpose(data.variables['nv'].data) - 1
            except AttributeError:
                self.trinodes = np.transpose(data.variables['nv'].data) - 1
        if 'triele' in datavar:
            try:
                setattr(self, 'triele', data.variables['triele'].data)
            except AttributeError:
                setattr(self, 'triele', data.variables['triele'])
        else:
            try:
                self.triele = np.transpose(data.variables['nbe'].data) - 1
            except AttributeError:
                self.triele = np.transpose(data.variables['nbe'].data) - 1

        # approximate bathymetry at the elements
        if self.threeD:
            if debug:
                print '\tapproximating bathymetry...'
            self.hele = (self.h[self.trinodes[0,:]-1] + \
                         self.h[self.trinodes[1,:]-1] \
                         + self.h[self.trinodes[2,:]-1]) / 3

        # load information relating to sigma layers and # of nodes/elements
        if self.threeD:                
            if debug:
                print '\tloading layers...'

            lsig, nodes = self.siglev.shape
            self.zlay = self.siglay[0:lsig-1, 0:1]
            self.zlev = self.siglev[0:lsig, 0:1]
            self.nlevel = self.siglay.shape[0]
            self.nsiglay = len(self.zlay)
            self.nsiglev = len(self.zlev)
        self.nele = self.lonc.shape[0]
        self.nnode = self.lon.shape[0]

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
        if debug:
            print '\tshifting grid to upper right cartesian...'
        self.x = self.x - self.xmin
        self.y = self.y - self.ymin

        # load time information, julian and matlab times
        try:
            self.mjd = data.variables['time'].data
        except AttributeError:
            # exception for nc.Dataset
            self.mjd = data.variables['time']

        self.gridtype = 'fvcom'
        
        if debug:
            print '...passed'


class load_scatter_grid:
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
            print 'retrieving data from ' + osp.basename(gridpath) + '...'
        if not osp.exists(gridpath):
            print '...the grid file was not found.'
            sys.exit()
        # enclose in try statement?
        data = sio.loadmat(gridpath)

        if debug:
            print 'loading grid variables...'
            
        # load grid variables
        # add __slots__?
        gridvars = ['lon', 'lat', 'time', 'uttc', 'uuss', 'uwnd', 'vttc', \
                    'vuss', 'vwnd']

        if debug:
            print '\tsetting attributes...'
        for key in gridvars:
            setattr(self, key, data[key])
        
        # define the lcc projection
        if debug:
            print '\tsetting lcc projection...'
            
        self._xmax = np.nanmax(self.lon)
        self._xmin = np.nanmin(self.lon)
        self._ymax = np.nanmax(self.lat)
        self._ymin = np.nanmax(self.lat)

        self._xavg = (self._xmax + self._xmin) * 0.5
        self._yavg = (self._ymax + self._ymin) * 0.5
        self._ylower = (self._ymax - self._ymin) * 0.25 + self._ymin
        self._yupper = (self._ymax - self._ymin) * 0.75 + self._ymin

        self._projstr = 'lcc +lon_0=' + str(self._xavg) + ' +lat_0=' \
                       + str(self._yavg) + ' +lat_1=' + str(self._ylower) \
                       + ' +lat_2=' + str(self._yupper)
        self.proj = pyproj.Proj(proj=self._projstr)

        self.mlon, self.mlat = np.meshgrid(self.lon, self.lat)
        self.lon, self.lat = self.mlon.flatten(), self.mlat.flatten()
        self.x, self.y = self.proj(self.lon, self.lat)
        self.mx, self.my = self.proj(self.mlon, self.mlat)

        # save the mask
        self.mask = np.isnan(self.uttc)[:,:,0]

        if debug:        
            print '\tflattening data fields...'

        # flatten the data as it is basically scattered in xy
        # self.uwndf = flattime(self.uwnd)
        # self.vwndf = flattime(self.vwnd)
        # self.uussf = flattime(self.uuss)
        # self.vussf = flattime(self.vuss)
            
        self.uttcf = flattime(self.uttc)
        self.vttcf = flattime(self.vttc)

        self.ustep1 = self.uttcf[0,:]
        self.vstep1 = self.vttcf[0,:]
        self.ustep2 = self.uttcf[0,:]
        self.vstep2 = self.vttcf[0,:]

        self.mjd = dn2mjd(self.time)
        self.gridtype = 'scatter'
        self.threeD = False
        
        if debug:
            print '...passed.'

        
class load_time_var:
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
            print 'loading time variables...'

        # Ensures internal step and output step are scalar multiples
        self.internalstep = np.float64(settings.lps)
        self.outputstep = np.float64(settings.ops)
        if not np.mod(self.internalstep, self.outputstep) == 0:
            print 'output and internal steps must be evenly divisible.'
            sys.exit()

        self.mjd = grid.mjd
        
        self.totalsteps = np.float64(settings.total)
        self.startstep = np.float64(settings.start)

        if debug:
            print '\tconverting measurements...'
        self.mdn = mjd2dn(self.mjd)
        self.ntime = self.mjd.shape[0]

        # convert mjd to jd, then to a datetime
        self.jd = mjd2jd(self.mjd)
        self._dates = jd2dt(self.jd)

        if debug:
            print '\tcalculating stepping...'
            
        # calculate stepping of input data, internal and output stepping
        self.instp = np.round((Time(self.jd[1], format='jd') \
                               -Time(self.jd[0], format='jd')).sec)

        self.dti = self.instp/self.internalstep
        self.dtout = self.instp/self.outputstep

        self.finishstep = self.startstep + self.totalsteps

        # add self.itout, iint, i2, int2, outt?

        # number of outputs to complete
        self.nouts = (self.totalsteps * self.outputstep) + 1

        if debug:
            print '...passed'


class load_part:
    """
    Loads the Lagrangian Particle subclass.

    --VARIABLES--
    nparts : number of particles
    locs : initial locations
    """
    def __init__(self, grid, ptime, settings, debug=False):
        if debug:
            print 'loading particle variables...'

        locs_file = settings.locs_file

        if debug:
            print '\tretrieving initial positions from ' \
                + osp.basename(locs_file) + '...'
        if not osp.exists(locs_file):
            print '...cannot find location file.'
            sys.exit()

        # set up initial positions, ignore missing data
        locs = np.genfromtxt(locs_file, comments='#', autostrip=True)
        
        if not len(settings.lon) == len(settings.lat):
            print 'number of longitudes and latitudes given must match.'
            sys.exit()

        # load lon/lat and elem given from settings
        if settings.lon and settings.lat:
            if debug:
                print '\tcollecting additional coordinates...'
                
            lon = np.asarray(settings.lon)
            lat = np.asarray(settings.lat)
            locs_opt = np.vstack((lon,lat)).T
            self.init_locs = np.vstack((locs,locs_opt))        

        else:
            self.init_locs = locs
            
        # initialize arrays for lagrangian particle(s)
        self.init_locs = self.init_locs[~np.isnan(self.init_locs).any(axis=1)]
        self.nparts = self.init_locs.shape[0]

        npts = self.nparts
        nouts = ptime.nouts
        out_path = settings.out_path
        # overwrite protection        
        if debug:
            print '\tfinding output file ' + osp.basename(out_path) + '...'
        if osp.isfile(out_path):
            choice = raw_input('\toutput file already exists. overwrite ' \
                               '(y/n)? ')
            if choice in ['yes', 'y', 'Y']:
                try:
                    os.remove(out_path)
                except OSError:
                    print 'could not overwrite at this time.'
                    sys.exit()
            elif choice in ['no', 'n', 'N']:
                print 'aborting...'
                sys.exit()
                
        # initialize netCDF object for writing
        try:
            if debug:
                print '\tcreating nc file...'
            ncid = nc.Dataset(out_path, 'w', format='NETCDF3_CLASSIC')
        except IOError:
            print 'could not write to output file at this time.'
            sys.exit()

        # create parameters
        if debug:
            print '\tsetting up dimensions and variables...'
        ncid.createDimension('time', nouts)
        ncid.createDimension('npts', npts)
        ncid.createVariable('x', 'd', ('time', 'npts'))
        ncid.createVariable('y', 'd', ('time', 'npts'))
        ncid.createVariable('lon', 'd', ('time', 'npts'))
        ncid.createVariable('lat', 'd', ('time', 'npts'))
        ncid.createVariable('u', 'd', ('time', 'npts'))
        ncid.createVariable('v', 'd', ('time', 'npts'))
        ncid.createVariable('time', 'd', ('time'))

        ncid.__setattr__('gridtype', grid.gridtype)
        ncid.__setattr__('history', 'created on ' + \
                         time.ctime(time.time()) + 'by PytoPlankton')

        if debug:
            print '\tcollecting initial step information...'
        # last step information
        self.lon0 = self.init_locs[:,0]
        self.lat0 = self.init_locs[:,1]
        # self.x0 = self.
        # self.y0 = self.

        # if self.threeD:
        #     self.h0 = 
        #     self.z0 =
        #     self.w0 =

        # create instances for immediate data (current step locations)
        self.xi = np.empty((npts,))
        self.yi = np.empty((npts,))
        self.loni = np.empty((npts,))
        self.lati = np.empty((npts,))
        self.ui = np.empty((npts,))
        self.vi = np.empty((npts,))

        if grid.threeD:
            self.hi = np.empty((npts,))
            self.zi = np.empty((npts,))
            self.wi = np.empty((npts,))

            ncid.createVariable('w', 'd', ('time','npts'))
            ncid.createVariable('z', 'd', ('time','npts'))
            ncid.createVariable('h', 'd', ('time','npts'))

        self.loop = 0

        ncid.variables['lon'][0,:] = self.init_locs[:,0]
        ncid.variables['lat'][0,:] = self.init_locs[:,1]
        # ncid.variables['x'][0,:] = self.x0
        # ncid.variables['y'][0,:] = self.y0
        # ncid.variables['u'][0,:]
        # ncid.variables['v'][0,:]

        # if self.threeD:
        #    ncid.variables['w'][0,:] =
        #    ncid.variables['z'][0,:] =
        #    ncid.variables['h'][0,:] =

        ncid.close()

        # difference between position / velocity and their absolutes?
        if debug:
            print '...passed'
