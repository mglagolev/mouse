#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from ordering_functions import *
import random

parser = argparse.ArgumentParser(description = 'Create a regular pattern of cubic domains')

parser.add_argument('--natom', type = int, nargs = 1, help = 'number of particles')

parser.add_argument('--size', type = float, nargs = 1, help = 'overall cell size')

parser.add_argument('--nbin', type = int, nargs = 1, help = 'number of bins in the histogram')

parser.add_argument('--rmin', type = float, nargs = 1, help = 'histogram start')

parser.add_argument('--rmax', type = float, nargs = 1, help = 'histogram end')

parser.add_argument('--icells', type = int, nargs = 3, help = 'number of cells: ix iy iz')

parser.add_argument('--pdb', type = str, nargs = '*', help = 'output .pdb file (optional)')

args = parser.parse_args()


def chequers(x, y, z, ix, iy, iz):
	nx = int( x * ix )
	ny = int( y * iy )
	nz = int( z * iz )
	if nx == 0 and ny == 0 and nz == 0:
		return random.choice(['A', 'B'])
	if nz % 2 == 0:
		if ny % 2 == 0:
			if nx % 2 == 0:
				return 'A'
			else:
				return 'B'
		else:
			if nx % 2 == 0:
				return 'B'
			else:
				return 'A'
	else:
		if ny % 2 == 0:
			if nx % 2 == 0:
				return 'B'
			else:
				return 'A'
		else:
			if nx % 2 == 0:
				return 'A'
			else:
				return 'B'

box = vector3d.Vector3d(args.size[0], args.size[0], args.size[0])
frame = Config()
frame.set_box(box)


for i in range(args.natom[0]):
	x = random.random()
	y = random.random()
	z = random.random()
	kind = chequers(x, y, z, args.icells[0], args.icells[1], args.icells[2])
	atom = Atom()
	pos = vector3d.Vector3d(x*args.size[0], y*args.size[0], z*args.size[0])
	atom.pos = pos
	atom.type = kind
	frame.insert_atom(atom)

if len(args.pdb) > 0:
	frame.write_pdb(args.pdb[0])

rdfs = CalculatePairCorrelationFunctions(frame, rmin = args.rmin[0], rmax = args.rmax[0], nbin = args.nbin[0])

pair_types = rdfs['data'].keys()[:]

pair_types.sort()

print rdfs['norm']

header = "# Distance        "

for pair_type in pair_types:
	header += pair_type.rjust(15) + "   "
print(header)

for i in range(len(rdfs['r'])):
	dataline = str(rdfs['r'][i]).rjust(15) + "   "
	for pair_type in pair_types:
		dataline += str(rdfs['data'][pair_type][i]/rdfs['v'][i]/rdfs['norm'][pair_type]).rjust(15) + "   "
	print(dataline)
