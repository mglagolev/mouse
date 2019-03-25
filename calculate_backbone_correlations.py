#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from backbone_correlation_functions import *

parser = argparse.ArgumentParser(description = 'Calculate static structure factor for one or more frames.')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--bondtypes', type = str, nargs = '+', help = 'types of bonds')

parser.add_argument('--max', type = int, nargs = '?', help = 'maximum distance between the bonds')

parser.add_argument('--chains', type = str, nargs = '+', help = 'chain numbers for analysis')

args = parser.parse_args()

for in_data in args.frames:
	print >> sys.stderr, '\r',
	frame = read_data_typeselect(in_data)
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+' Frame nbonds: '+str(frame.n_bond())+'\n')
	if len(args.chains) > 0:
		frame = select_chains(frame, args.chains)
		CalculateBackboneCorrelations(frame, args.bondtypes, args.max)
	else:
		CalculateBackboneCorrelations(frame, args.bondtypes, args.max)
print '\n',
