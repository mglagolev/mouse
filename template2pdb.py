#!/usr/bin/env python

import sys
import os
import argparse
from lammpack_types import Config
from pla_reverse import *

parser = argparse.ArgumentParser(description = 'Convert reverse mapping template to .pdb format')

parser.add_argument('template_file', metavar = '.dat', type = str, nargs = 1, help = 'Reverse mapping template, dict format')

parser.add_argument('pdb_file', metavar = '.pdb', type = str, nargs = 1, help = 'Template in PDB format')

args = parser.parse_args()

config = templateToConfig(args.template_file[0])
config.write_pdb(args.pdb_file[0])
