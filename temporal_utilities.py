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
import seaborn as sns
from astropy.time import Time

# local imports


def dn2dt(dn):
    """
    Converts a Matlab datenum to a Python datetime.
    """
    return datetime.fromordinal(int(dn)) + timedelta(days=dn%1) \
            - timedelta(days = 366)

    
def jd2dt(jd):
    """
    Converts a Julian Time to a Python datetime.
    """
    t = Time(jd, scale='utc', format='jd')
    return t.datetime


def mjd2dn(mjd):
    """
    Converts a Modified Julian Time to a Matlab datenum.
    """
    return mjd[:] + 678942.0


def mjd2jd(mjd):
    """
    Converts a Modified Julian Time to a Julian Time.
    """
    return mjd[:] + 2400000.5

    
def jd2isot(jd):
    """
    Converts a Julian Time to an ISO 8601 compliant timestamp.
    """
    t = Time(jd, scale='utc', format='jd')
    return t.isot
