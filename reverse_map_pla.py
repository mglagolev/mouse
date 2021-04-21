#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from pla_reverse import *

parser = argparse.ArgumentParser(description = 'Create atomistic .gro file from coarse-grained .gro file for A-graft-B model of PLA')

parser.add_argument('in_gro', metavar = '.GRO', type = str, nargs = 1, help = 'GROMACS coordinate file')

parser.add_argument('--structure', metavar = 'STRUCTURE_FILE', type = str, nargs = '+', help = 'dictionary files with structure of start, middle and end unit of PLA')

parser.add_argument('--out', metavar = '.GRO', type = str, nargs = 1, help = 'GROMACS coordinate file')

parser.add_argument('--scaling', metavar = "SCALING_FACTOR", type = float, nargs = '+', help = "Use internally calculated scaling factor with this weight; s = (s0 - 1) * SCALING_FACTOR + 1. In case two values are provided the first one is used for downscaling, the second for upscaling", default = 1.0)

parser.add_argument('--adjust-pendant', help = "Adjust the position of the methyl group according to positon of the coarse-grained B pendant", action = 'store_true')

parser.add_argument('--side-scaling', metavar = "PENDANT SCALING FACTOR", type = float, nargs = '+', help = "Use pendant scaling factor with this weight; s = (s0 - 1) * PENDANT SCALING_FACTOR + 1. In case two values are provided the first one is used for downscaling, the second for upscaling", default = 1.0)


args = parser.parse_args()

if len(args.scaling) == 1: scaling = [args.scaling, args.scaling]
elif len(args.scaling) == 2: scaling = [args.scaling[0], args.scaling[1]]

if len(args.side_scaling) == 1: side_scaling = [args.side_scaling, args.side_scaling]
elif len(args.side_scaling) == 2: side_scaling = [args.side_scaling[0], args.side_scaling[1]]

print "ADJUST PENDANT ", args.adjust_pendant
for in_data in args.in_gro:
	print >> sys.stderr, '\r',
	frame = read_data_typeselect(in_data)
	sys.stderr.write('Frame natom: '+str(frame.n_atom())+'\n')
	aa_config = pla_reverse_map(frame, args.structure, scaling, side_scaling, args.adjust_pendant)
	aa_config.write_gro(args.out[0])
