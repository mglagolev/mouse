#!/usr/bin/env python3

import sys
import os
import argparse
from lammpack_types import Config
from ordering_functions import *

parser = argparse.ArgumentParser(description = 'Calculate local values of orientational order parameter based on backbone bonds.')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--bondtypes', type = str, nargs = '+', help = 'types of bonds')

parser.add_argument('--max', type = float, nargs = '?', help = 'maximum distance between the bonds')

parser.add_argument('--min', type = float, nargs = '?', help = 'minimum distance between the bonds')

parser.add_argument('--chains', type = str, nargs = '+', help = 'chain numbers for analysis')

parser.add_argument('--mode', type = str, nargs = '?', help = 'output mode: average | histo', default = 'average')

parser.add_argument('--histo-min', type = float, nargs = '?', help = 'histogram minimum', default = -0.5)

parser.add_argument('--histo-max', type = float, nargs = '?', help = 'histogram maximum', default = 1.)

parser.add_argument('--histo-bins', type = int, nargs = '?', help = 'number of bins', default = 75)

parser.add_argument('--subcell', type = float, nargs = '+', help = 'rectangular subcell boundaries, normalized')

parser.add_argument('--reference-residue', type = str, nargs = '+', help = 'residue types of atoms for bonds used as reference')

parser.add_argument('--same-molecule', action = 'store_true', help = 'take into accounts bond in the same molecule')

parser.add_argument('--out-pdb', type = str, nargs = '*', help = "output .pdb with atom types encoding local ordering value, e.g. CXX, where XX is crystallinity value in %")

args = parser.parse_args()

for in_data in args.frames:
	sys.stderr.write("\r")
	if not args.same_molecule: readOptions = {'pdb' : { 'assignMolecules' : {'type' : 'resSeq', 'cluster' : True } } }
	else: readOptions = {}
	frame = read_data_typeselect(in_data, options = readOptions)
	sys.stderr.write('Read frame (natom): '+str(frame.n_atom())+'\n')
	try:
		if len(args.chains) > 0:
			frame = select_chains(frame, args.chains)
	except TypeError:
		pass
	sys.stderr.write('After molecules selection: '+str(frame.n_atom())+'\n')
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
		outPdb = args.out_pdb[0]
		storeAsAtomtypes = True
	except: storeAsAtomtypes = False
	if len(args.subcell) > 0:
		frame = select_rectangular(frame, cutx, xrelmin, xrelmax, cuty, yrelmin, yrelmax, cutz, zrelmin, zrelmax)
	sys.stderr.write('After region selection: '+str(frame.n_atom())+'\n')
	s = CalculateOrientationOrderParameter(frame, args.bondtypes, args.min, args.max, mode = args.mode, storeAsAtomtypes = storeAsAtomtypes, smin = args.histo_min, smax = args.histo_max, sbins = args.histo_bins, referenceResTypes = args.reference_residue, sameMolecule = args.same_molecule)
	if storeAsAtomtypes:
		frame.write_pdb(outPdb, hide_pbc_bonds = True)
	if args.mode == 'average': sys.stdout.write(str(s) + "\n")
	elif args.mode == 'histo': sys.stdout.write(printHistogram(s, norm = False))
		
