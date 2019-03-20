#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from spatial_functions import *

parser = argparse.ArgumentParser(description = 'Calculate average helix tube radius of the macromolecules.')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--atomtypes', type = str, nargs = '+', default=["All"], help = 'types of atoms considered for calculation')

parser.add_argument('--chains', type = int, nargs = '+', help = 'chain numbers for analysis')

parser.add_argument('--period', type = float, nargs = 1, help = 'number of units per helix turn')

args = parser.parse_args()

for in_data in args.frames:
	print >> sys.stderr, '\r',
	frame = read_data_typeselect(in_data)
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+'\n')
	if len(args.chains) > 0:
		frame = select_chains(frame, args.chains)
	frame.convertCoordsCutToUncut()
	average_tube_radius = CalculateHelixTubeRadius(frame, args.atomtypes, args.period[0])
	print average_tube_radius
