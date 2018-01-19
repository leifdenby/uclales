# -*- coding: utf-8 -*-
"""
Create density weighted uniform particle distribution for the Lagrangian
Particle Tracking Module in UCLA-LES. As UCLA-LES uses the anelastic    
equations, we need to start with a density weighted distribution.       

Leif Denby 19/2/2018, Based on script by Bart van Stratum Aug 2012
"""
from __future__ import print_function

from math import floor
import argparse
import os

import namelist_python

class Defaults:
    tstart      = 0.          # release time of particles
    nxy         = 128**2.     # Goal number of particles per level
    output_filename = "partstartpos"


argparser = argparse.ArgumentParser(
    __doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
argparser.add_argument(
    'namelist', type=str,
    help="UCLALES namelist file"
)
argparser.add_argument(
    '--tstart', type=float, default=Defaults.tstart,
    help="release time of particles"
)
argparser.add_argument(
    '--nxy', type=float, default=Defaults.nxy,
    help="number of particles per level"
)

args = argparser.parse_args()

namelist_filename = args.namelist
tstart = args.tstart
nxy = args.nxy

if not os.path.exists(namelist_filename):
    raise Exception("Couldn't find UCLALES namelist file `{}`".format(
        namelist_filename
    ))
else:
    nml = namelist_python.read_namelist_file(namelist_filename)

    th00 = nml.data.model.th00
    ps = nml.data.model.ps[0]

    dx = nml.data.model.deltax
    dy = nml.data.model.deltay
    dz = nml.data.model.deltaz

    # need to exclude ghost cells
    nx = nml.data.model.nxp - 2
    ny = nml.data.model.nyp - 2
    nz = nml.data.model.nzp - 2

    if nx != ny or dx != dy:
        raise NotImplementedError("Cannot work with non-uniform xy-grid")

    xysize = nx*dx

    dzrat = nml.data.model.dzrat
    if dzrat != 1.0:
        raise NotImplementedError("Cannot work with stretched grid")

    p00 = 1.0e5  # reference pressure, hard-coded in UCLALES



# No need to change below
# ----------------------------------------------
# CONSTANTS:
cp          = 1005.
R           = 287.04
g           = 9.8
p00i        = 1. / p00
rcp         = R/cp
cpr         = cp/R

# Calculate number of particles.
# Required before writing the output, first line should
# contain the total number of particles.
npart       = 0
pi00        = cp * (ps * p00i)**rcp + g * (0.5 * dz) / th00
for k in range(nz):
  dz0       = -float(k+1)  * dz
  pi0       = pi00 + g * dz0 / th00
  dn0       = ((cp**(1.-cpr)) * p00) / (R * th00 * pi0**(1.-cpr))
  nxyl      = round((nxy * dn0)**0.5)
  npart    += nxyl**2.

# Write to partstartpos
with open(Defaults.output_filename,'w') as partfile:
  partfile.write(str(int(npart))+'\n')
  z = dz / 2.
  for k in range(nz):
    dz0       = -float(k+1)  * dz
    pi0       = pi00 + g * dz0 / th00
    dn0       = ((cp**(1.-cpr)) * p00) / (R * th00 * pi0**(1.-cpr))
    nxyl      = round((nxy * dn0)**0.5)
    dxy       = xysize / float(nxyl)
    x         = dxy / 2.
    y         = dxy / 2.
    print('k = %3i / %3i, z = %.2f, # = %.3f' % (k+1,nz,z,nxyl), end='\r')
    #print '%3i  %.2f  %.3f' % (k+1,z,nxyl)
    for i in range(int(nxyl)):
      for j in range(int(nxyl)):
        partfile.write('%i  %.3e  %.3e  %.3e \n'%(tstart,x,y,z))
        y += dxy
      y = dxy/2.
      x += dxy
    z += dz

print("Particle start positions written to {fn}".format(
    fn=Defaults.output_filename)
)
