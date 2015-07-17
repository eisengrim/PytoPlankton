v#!/usr/bin/env python
# encoding: utf-8

# library imports
from __future__ import division
import sys
import os
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

# local imports
from temporal_utilities import *
from spatial_utilities import *

M_STAGE = 4
A = [0, 0.5, 0.5, 1]
B = [1/6, 1/3, 1/3, 1/6]
C = [0, 0.5, 0.5, 1]


def solver(particle, grid, time, settings):
    """
    Calls the numerical solver to be used for moving the Lagrangian
    particle(s).
    """
    if settings['solver'] == 'rk4':
        _rungekutta4(particle, grid, time, settings)
    elif settings['solver'] == 'fwd_euler':
        _forward_euler(particle, grid, time, settings)
    else:
        print '{} is not a valid solver'.format(settings['solver']
        sys.exit()


def _rungekutta(particle, grid, time, settings):
    """
    Executes an RK4 numerical scheme. Integrates all particle positions using
    the velocity fields given from the grid.
    """

    chi_x = np.zeros((particle.nparts,4))
    chi_y = np.zeros((particle.nparts,4))

    for i in xrange(0, M_STAGE):
        xpt  = lag['xl']  + (A[i]*lag['timestep'])*chix[:,i]
        ypt  = lag['yl']  + (A[i]*lag['timestep'])*chiy[:,i]

        uin  = ((1-c_rk[ns])*grid['ustep1'] + c_rk[ns]*grid['ustep2'])
        vin  = ((1-c_rk[ns])*grid['vstep1'] + c_rk[ns]*grid['ustep2'])

        usam=interp(settings,grid,uin,xpt,ypt)
        vsam=interp(settings,grid,vin,xpt,ypt)

        chix[:,ns] = usam
        chiy[:,ns] = vsam
 
    xpt=lag['xl']
    ypt=lag['yl']
    for ns in range(0,mstage):
        xpt = xpt + lag['timestep']*b_rk[ns]*chix[:,ns]
        ypt = ypt + lag['timestep']*b_rk[ns]*chiy[:,ns]

    lag['xn']=xpt
    lag['yn']=ypt
    lag['u']=interp(settings,grid,uin,lag['xn'],lag['yn'])
    lag['v']=interp(settings,grid,vin,lag['xn'],lag['yn'])

    return


def _forward_euler(particle, lag, time, settings):
    """
    Executes a Forward Euler method, integrating all particle positions using
    the fields given from the grid.
    """
    pass
