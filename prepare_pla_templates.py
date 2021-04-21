#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from pla_reverse import *

parser = argparse.ArgumentParser(description = 'Create templates for reverse mapping of A-graft-B model of PLA from atomistic .gro file')

parser.add_argument('in_gro', metavar = '.GRO', type = str, nargs = 1, help = 'GROMACS coordinate file')

parser.add_argument('--types', metavar = 'UNIT_TYPE', type = str, nargs = '+', help = 'unit types (SLS, SLM, SLE)')

parser.add_argument('--name', metavar = 'TEMPLATE_NAME', type = str, nargs = 1, help = 'base name for the templates')

args = parser.parse_args()

for unitType in args.types:
	print >> sys.stderr, '\r',
	frame = read_data_typeselect(args.in_gro[0])
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+'\n')
	writeUnitStructures(frame, unitType, args.name[0])
