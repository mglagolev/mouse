#!/usr/bin/env python3

import sys
import os
import argparse
from lammpack_types import Config
from ordering_functions import *

parser = argparse.ArgumentParser(description = 'Calculate radial distribution functions')

parser.add_argument('frame', metavar = 'data', type = str, nargs = 1, help = 'snapshot file')

parser.add_argument('--rmin', type = float, nargs = '?', default = 0., help = 'minimum distance')

parser.add_argument('--rmax', type = float, nargs = '?', default = 0., help = 'maximum distance')

parser.add_argument('--nbin', type = int, nargs = '?', default = 100, help = 'number of bins')

parser.add_argument('--no-norm', action = 'store_false', help = 'do not normalize the distributions')

args = parser.parse_args()

sys.stderr.write('\r')
frame = read_data_typeselect(args.frame[0])
sys.stderr.write('Frame natom: '+str(frame.n_atom())+' Frame nbonds: '+str(frame.n_bond())+'\n')

rdfs = CalculateComRdfs(frame, args.rmin, args.rmax, args.nbin)

printRdfs(rdfs, args.no_norm)

