#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from lammpack_misc import *
from clustering_functions import *

parser = argparse.ArgumentParser(description = 'Calculate atom lists for every molecule in the configuration')

parser.add_argument('frames', metavar = 'frame', type = str, nargs = '+', help = 'snapshot file')

args = parser.parse_args()

for in_data in args.frames:
	print >> sys.stderr, '\r',
	frame = read_data_typeselect(in_data)
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+' Frame nbonds: '+str(frame.n_bond())+'\n')
	makeNeighborlistsFromBonds(frame)
	aggregates = determineClusters(frame)
	for aggregate in aggregates:
		for atom in aggregate:
			print(atom.type + " " + str(atom.num) + "   "),
		print "\n\n"
print '\n',
