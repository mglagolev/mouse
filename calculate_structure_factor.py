#!/usr/bin/env python3

import sys
import os
import argparse
from lammpack_types import Config
from lammpack_misc import read_data_typeselect, select_chains, selectByAtomFields
import structure_factor
import math

parser = argparse.ArgumentParser(description = 'Calculate structure factor')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--g', type = int, nargs = '+', help = 'number of points', default = 50)

parser.add_argument('--ksqmax', type = int, nargs = '+', help = 'maximum k^2 value', default = 2500)

parser.add_argument('--atomtypes', type = str, nargs = '+', help = 'types of atoms')

parser.add_argument('--chains', type = str, nargs = '+', help = 'chain numbers for analysis')

args = parser.parse_args()

for in_data in args.frames:
	sys.stderr.write("\r")
	frame = read_data_typeselect(in_data, options = {"lammps_data":{"bonds":False, "angles":False, "dihedrals":False}})
	if len(args.atomtypes) > 0:
		frame = selectByAtomFields(frame, { "type" : args.atomtypes })
		sys.stderr.write('Selected by type: '+str(frame.n_atom())+'\n')
	sys.stderr.write('Read frame (natom): '+str(frame.n_atom())+'\n')
	try:
		if len(args.chains) > 0:
			frame = select_chains(frame, args.chains)
	except TypeError:
		pass
	sys.stderr.write('After molecules selection: '+str(frame.n_atom())+'\n')
	sk, ns = structure_factor.calculateStructureFactor(frame, args.g[0], args.ksqmax[0])
	sys.stdout.write("#k		sk/ns\n")
	for i in range(len(sk)):
	    if ns[i] > 0:
	        sys.stdout.write(str(math.sqrt(i) / min(frame.box().x, frame.box().y, frame.box().z) ) 	+ "	" + str(sk[i]/float(ns[i])) + "\n")
	    
