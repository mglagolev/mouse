#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from ssf_fft import skfft, skfinalize

parser = argparse.ArgumentParser(description = 'Calculate static structure factor for one or more frames.')

parser.add_argument('frames', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--data-type', type = str, nargs = 1, help = 'input data type: "lammps_data" or "pdb"')

parser.add_argument('--out', type = str, nargs = 1, help = 'output file')


args = parser.parse_args()

for in_data in args.frames:
	print >> sys.stderr, '\r',
	frame = Config()
	frame.read_lmp_data(in_data)
	skfft(frame)
	print >> sys.stderr, 'Processed', os.path.basename(in_data),
	
skfinalize(frame, args.out[0])

print '\n',
