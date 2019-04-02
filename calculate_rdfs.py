#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from ordering_functions import *

def printRdfs(rdfs):
	nlines = len(rdfs['r'])
	sys.stdout.write("#r	vol")
	orderedKeys = []
	for key in rdfs['data']:
		orderedKeys.append(key)
	orderedKeys.sort()
	for key in orderedKeys:
		sys.stdout.write("	"+str(key))
	sys.stdout.write("\n")
	for nline in range(nlines):
		r, v = rdfs['r'][nline], rdfs['v'][nline]
		sys.stdout.write(str(r))
		sys.stdout.write("	" + str(v))
		for key in orderedKeys:
			sys.stdout.write("	" + str(rdfs['data'][key][nline] / rdfs['norm'][key] / r / v ))
		sys.stdout.write("\n")

parser = argparse.ArgumentParser(description = 'Calculate radial distribution functions')

parser.add_argument('frame', metavar = 'data', type = str, nargs = 1, help = 'snapshot file')

parser.add_argument('--atomtypes', type = str, nargs = '+', help = 'types of atoms')

parser.add_argument('--rmin', type = float, nargs = '?', default = 0., help = 'minimum distance')

parser.add_argument('--rmax', type = float, nargs = '?', default = 0., help = 'maximum distance')

parser.add_argument('--nbin', type = int, nargs = '?', default = 100, help = 'number of bins')

parser.add_argument('--chains', type = str, nargs = '+', help = 'chain numbers for analysis')

parser.add_argument('--same-molecule', action = 'store_true', help = 'calculate for atoms in the same molecule')

args = parser.parse_args()

print >> sys.stderr, '\r',
frame = read_data_typeselect(args.frame[0])
sys.stderr.write('Frame natom: '+str(frame.n_atom())+' Frame nbonds: '+str(frame.n_bond())+'\n')
try:
	if len(args.chains) > 0:
		frame = select_chains(frame, args.chains)
		rdfs = CalculatePairCorrelationFunctions(frame, args.atomtypes, args.rmin, args.rmax, args.nbin, args.same_molecule)
except TypeError:
	rdfs = CalculatePairCorrelationFunctions(frame, args.atomtypes, args.rmin, args.rmax, args.nbin, args.same_molecule)

printRdfs(rdfs)


