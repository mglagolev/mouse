#!/usr/bin/env python3

import sys
from lammpack_types import Config
import argparse
import vector3d

parser = argparse.ArgumentParser(description = 'Convert LAMMPS data file into PDB')
parser.add_argument('input', action = "store", help = "input file")
parser.add_argument('output', action = "store", help = "output file")
parser.add_argument("--center", action = "store_true", help = "Center the box around zero")
args = parser.parse_args()

frame = Config()

frame.read_lmp_data(args.input)

if args.center:
	frame.center_box()

frame.write_pdb(args.output)
