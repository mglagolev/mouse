#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from ordering_functions import *

parser = argparse.ArgumentParser(description = 'Calculate average value of crystallinity parameter (3cos^2(theta)-1)/2.')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--bondtypes', type = str, nargs = '+', help = 'types of bonds')

parser.add_argument('--chains', type = str, nargs = '*', help = 'chain numbers for analysis')

parser.add_argument('--subcell', type = float, nargs = '*', help = 'rectangular subcell boundaries, normalized xmin, xmax, [ymin, ymax, [zmin, zmax]]]')

parser.add_argument('--reference', type = float, nargs = '*', help = 'optional reference direction vector components')

parser.add_argument('--mol-id', type = str, nargs = 1, default = 'resSeq', help = ".pdb value to identify individual molecules, resSeq (default), or chainId")

args = parser.parse_args()

readOptions = { "pdb" : { "assignMolecules" : { "type" : args.mol_id, "cluster" : True }}}

for in_data in args.frames:
	print >> sys.stderr, '\r',
	frame = read_data_typeselect(in_data, options = readOptions)
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+'\n')
	try:
		if len(args.chains) > 0:
			frame = select_chains(frame, args.chains)
	except TypeError:
		pass
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+'\n')
	xrelmin, xrelmax, yrelmin, yrelmax, zrelmin, zrelmax = 0., 1., 0., 1., 0., 1.
	cutx, cuty, cutz = False, False, False
	try:
		xrelmin = args.subcell[0]
		xrelmax = args.subcell[1]
		cutx = True
		yrelmin = args.subcell[2]
		yrelmax = args.subcell[3]
		cuty = True
		zrelmin = args.subcell[4]
		zrelmax = args.subcell[5]
		cutz = True
	except:
		pass
	try:
		refx = args.reference[0]
		refy = args.reference[1]
		refz = args.reference[2]
		reference_vector = vector3d.Vector3d(refx, refy, refz)
	except:
		reference_vector = vector3d.Vector3d(0., 0., 0.)
	if cutx or cuty or cutz:
		frame = select_rectangular(frame, cutx, xrelmin, xrelmax, cuty, yrelmin, yrelmax, cutz, zrelmin, zrelmax)
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+'\n')
	s = CalculateCrystallinityParameter(frame, args.bondtypes, reference_vector = reference_vector)
	print s
