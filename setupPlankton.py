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


class _load_settings:
    """
    Loads the settings for the Plankton object, and stores it in a subclass.
    It is necessary to change these parameters before running the Lagrangian
    tracker code. This is done by simply calling the Settings subclass and
    assigning a new value to the key.
    """
    def __init__(self, **kwargs):
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

    
class _load_grid(raw_data, settings, debug=False):
    """
    Loads the grid data for the subclass Grid.
    """
    pass


class _load_time_var(raw_data, settings, debug=False):
    """
    Loads the time variables subclass Time.
    """
    pass


class _load_part(raw_data, settings, debug=False):
    """
    Loads the Lagrangian Particle subclass.
    """
    pass
