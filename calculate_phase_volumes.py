#!/usr/bin/env python3

import sys
from lammpack_types import Config
import argparse
import vector3d
import spatial_functions

parser = argparse.ArgumentParser(description = 'Calculate phase volume in a system')
parser.add_argument('input', action = "store", help = "input file")
parser.add_argument('--ntrials', type = int, nargs = '?', default = 8000, help = 'number of trials')
parser.add_argument('--rprobe', type = float, nargs = '?', default = -1, help = 'probe radius')
args = parser.parse_args()

frame = Config()

frame.read_lmp_data(args.input)

phases, stdev = spatial_functions.CalculatePhaseVolume(frame, args.ntrials, args.rprobe)

for atype in sorted(phases):
	sys.stdout.write(str(atype) + "	" + str(phases[atype]) + "	" + str(stdev[atype]) + "\n")
