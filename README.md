# PytoPlankton v0.1

PytoPlankton is a particle-tracking software that uses Lagrangian modelling for
output from FVCOM models (as netCDF4) and MATLAB scatter data. The program is a 
Python implementation
of MATLAB code written by Mitchell O'Flaherty (BIO), based on a FORTRAN program
originally designed to simulate the motion of sea lice (the name PytoPlankton was 
inspired from the planktonic behaviour that sea lice exhibit while searching for 
hosts to attach to). Parts of PytoPlankton are adapted from Thomas Roc's PySeidon,
a validation suite for FVCOM outputs.Several scientific Python modules are called 
in PytoPlankton. These can be downloaded via Anaconda, a free Python distribution:

http://continuum.io/downloads

Both PySeidon and the original MATLAB code, LagTracker, can be found here:

https://github.com/GrumpyNounours/PySeidon

https://github.com/moflaher/lagtracker

