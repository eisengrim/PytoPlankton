v#!/usr/bin/env python
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
from temporal_utilities import *
from spatial_utilities import *


A_RK = [0, 0.5, 0.5, 1]
B_RK = [1/6, 1/3, 1/3, 1/6]
C_RK = [0, 0.5, 0.5, 1]


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
        sys.exit('{} is not a valid solver'.format(settings['solver'])


def _rungekutta(particle, grid, time, settings):
    """
    Executes an RK4 numerical scheme. Integrates all particle positions using
    the velocity fields given from the grid.
    """
    pass


def _forward_euler(particle, lag, time, settings):
    """
    Executes a Forward Euler method, integrating all particle positions using
    the fields given from the grid.
    """
    pass
