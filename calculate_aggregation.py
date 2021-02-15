#!/usr/bin/env python3

import sys
import os
import argparse
from lammpack_types import Config
from ordering_functions import *

parser = argparse.ArgumentParser(description = 'Determine the aggregates in the systems and print out the atoms belonging to each of them.')

parser.add_argument('frames', metavar = 'ATOMIC COORDINATES', type = str, nargs = '+', help = 'data file')

parser.add_argument('--atomtypes', type = str, nargs = '+', help = 'types of atoms')

parser.add_argument('--threshold', type = float, nargs = '?', default = 1.2, help = 'maximum distance between the atoms to be considered as clustered')

args = parser.parse_args()

for in_data in args.frames:
	sys.stderr.write('\r')
	frame = read_data_typeselect(in_data)
	sys.stderr.write('Read frame (natom): '+str(frame.n_atom())+'\n')
	frame = selectByAtomFields(frame, { "type" : args.atomtypes })
	sys.stderr.write('Selected by type: '+str(frame.n_atom())+'\n')
	frame.reassign_inums()
	sys.stderr.write('Box size: x ' + str(frame.box().x) + ' y ' + str(frame.box().y) + ' z ' + str(frame.box().z) + '\n')
	makeNeighborlistsFromDistances(frame, threshold = args.threshold)
	all_aggregates = determineClusters(frame)
	for aggregate in all_aggregates:
		atom_nums_array = []
		for atom in aggregate:
			atom_nums_array.append(atom.num)
		atom_nums_array.sort()
		sys.stdout.write(str(atom_nums_array) + "\n")
		
