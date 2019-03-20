#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from spatial_functions import *

parser = argparse.ArgumentParser(description = 'Calculate average end to end length of the macromolecules.')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--atomtypes', type = str, nargs = '+', default=["All"], help = 'types of atoms considered for calculation')

parser.add_argument('--chains', type = int, nargs = '+', help = 'chain numbers for analysis')

args = parser.parse_args()

for in_data in args.frames:
	print >> sys.stderr, '\r',
	frame = read_data_typeselect(in_data)
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+'\n')
	if len(args.chains) > 0:
		frame = select_chains(frame, args.chains)
	average_end2end = CalculateEnd2EndDistance(frame, args.atomtypes)
	print average_end2end
