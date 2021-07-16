#!/usr/bin/env python3

import sys
import os
import argparse
from lammpack_types import Config
from lammpack_misc import *
import end_to_end

parser = argparse.ArgumentParser(description = 'Calculate end-to-end distance of the blocks')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--atomtypes', type = str, nargs = '+', help = 'types of atoms')

parser.add_argument('--min', type = float, nargs = '?', help = 'histogram min', default = 0.)

parser.add_argument('--max', type = float, nargs = '?', help = 'histogram max', default = 30.)

parser.add_argument('--nbins', type = int, nargs = '?', help = 'histogram nbins', default = 300)

parser.add_argument('--chains', type = str, nargs = '+', help = 'chain numbers for analysis')

args = parser.parse_args()

for in_data in args.frames:
	sys.stderr.write("\r")
	frame = read_data_typeselect(in_data, options = {"lammps_data":{"bonds":False, "angles":False, "dihedrals":False}})
	sys.stderr.write('Read frame (natom): '+str(frame.n_atom())+'\n')
	try:
		if len(args.chains) > 0:
			frame = select_chains(frame, args.chains)
	except TypeError:
		pass
	sys.stderr.write('After molecules selection: '+str(frame.n_atom())+'\n')
	if args.atomtypes is not None:
		frame = selectByAtomFields(frame, { "type" : args.atomtypes })
	histograms = {}
	for atomtype in set(atom.type for atom in frame.atoms()):
		current_config = selectByAtomFields(frame, { "type" : [atomtype]})
		current_config = frame
		e2e_list = end_to_end.end_to_end_by_molecule(current_config)
		histograms[atomtype] = createHistogram(args.min, args.max, args.nbins)
		for e2e in e2e_list:
			updateHistogram(e2e, histograms[atomtype])
		normHistogram(histograms[atomtype])
		sys.stdout.write(printHistogram(histograms[atomtype]))
		
		
