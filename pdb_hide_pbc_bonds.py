#!/usr/bin/env python

import sys
from lammpack_types import Config


frame = Config()

frame.read_pdb(sys.argv[1])

frame.write_pdb(sys.argv[2], True)
