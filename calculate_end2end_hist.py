#!/usr/bin/env python3

import sys
import os
import argparse
from lammpack_types import Config
from end_to_end import *

parser = argparse.ArgumentParser(description = 'Calculate end-to-end distance of the blocks')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--type1', type = str, nargs = '+', help = 'atom type for block 1', default = "1")

parser.add_argument('--type2', type = str, nargs = '+', help = 'atom type for block 2', default = "2")

parser.add_argument('--chains', type = str, nargs = '+', help = 'chain numbers for analysis')

args = parser.parse_args()

for in_data in args.frames:
	sys.stderr.write("\r")
	frame = read_data_typeselect(in_data)
	sys.stderr.write('Read frame (natom): '+str(frame.n_atom())+'\n')
	try:
		if len(args.chains) > 0:
			frame = select_chains(frame, args.chains)
	except TypeError:
		pass
	sys.stderr.write('After molecules selection: '+str(frame.n_atom())+'\n')
	EndToEndHist(frame)
	sys.stdout.write("Result: " + "\n")
		
