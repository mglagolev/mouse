#!/usr/bin/python
import sys
import numpy as np
from numpy import *
factor = float(sys.argv[2])
s = open(sys.argv[1], 'r').read()
coords = eval(s)
for key in coords:
	coords[key] *= factor
o = open(sys.argv[3], 'w')
o.write(str(coords))
o.close()
