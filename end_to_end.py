#!/usr/bin/env python3

from lammpack_misc import *
from lammpack_types import *
from ordering_functions import createNumpyBondsArrayFromConfig
from ordering_functions import createNumpyAtomsArrayFromConfig
import vector3d
import numpy as np
from numpy import linalg as LA
import sys
import pdb
from lammpack_misc import GetMoleculeAtomlists
         
def end_to_end(atoms, sorting = "num"):
    """ Order the atoms based on the sorting attribute, and calculate the vector from firt to last atom """
    if len(atoms) == 0: raise NameError("Can not calculate end-to-end value for an empty list")
    atoms.sort(key = lambda x: getattr(x, sorting))
    return atoms[-1].pos - atoms[0].pos
    
def end_to_end_by_molecule(config, sorting = "num"):
    """ Calculate PBC-trimmed vectors from the first to last atom belonging to the same molecule """
    e2e_distances = []
    _, atomlists = GetMoleculeAtoms(config)
    for atomlist in atomlists:
        e2e_vector = vector_pbc_trim(end_to_end(atomlist, sorting = sorting), config.box())
        e2e_distances.append(e2e_vector.length())
    return e2e_distances
        
