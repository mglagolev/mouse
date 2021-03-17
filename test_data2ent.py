#!/usr/bin/env python3

import sys
from lammpack_types import Config
import argparse
import vector3d
import lammpack_misc

parser = argparse.ArgumentParser(description = 'Convert LAMMPS data file into PDB')
parser.add_argument('input', action = "store", help = "input file")
parser.add_argument('output', action = "store", help = "output file")
parser.add_argument("--center", action = "store_true", help = "Center the box around zero")
parser.add_argument("--hide-pbc-bonds", action = "store_true", help = "Hide the bonds transversing the periodic boundary conditions")
parser.add_argument('--atomtypes', type = str, nargs = '+', help = 'types of atoms')

args = parser.parse_args()

frame = Config()

frame.read_lmp_data(args.input, angles = False, dihedrals = False)

if args.atomtypes is not None:
	frame = lammpack_misc.selectByAtomFields(frame, { "type" : args.atomtypes })

if args.center:
	frame.center_box()

frame.write_pdb(args.output, args.hide_pbc_bonds)
