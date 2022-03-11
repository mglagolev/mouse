#!/usr/bin/env python3

import sys
import os
import argparse
from lammpack_types import *
from lammpack_misc import read_data_typeselect

parser = argparse.ArgumentParser(description = 'Enumerate the blocks of a multiblock copolymer and assign different mol_ids')

parser.add_argument('frame', metavar = 'data', type = str, nargs = 1, help = 'snapshot file')

parser.add_argument('output', metavar = 'data', type = str, nargs = 1, help = 'output')

parser.add_argument('--segment', metavar = 'number', type = int, nargs = 1, help = 'segment length for analysis')

args = parser.parse_args()

sys.stderr.write('\r')
frame = read_data_typeselect(args.frame[0])

atoms = frame.atoms()
atoms.sort(key=lambda atom: atom.num)
prevtype = "None"
mol_id = 0
for atom in atoms:
    atomtype = atom.type
    if atomtype != prevtype:
        mol_id += 1
    atom.mol_id = str(mol_id)
    prevtype = atom.type

n_bond_types = len(list(set(bond.type for bond in frame.bonds())))    
bond_num = len(frame.bonds())
for i in range(1, len(frame.atoms()) - args.segment[0] + 2):
    atom1 = frame.atom_by_num(i)
    atom2 = frame.atom_by_num(i + args.segment[0] - 1)
    if atom1.mol_id == atom2.mol_id:
        atom_type = atom1.type
        bond = Bond()
        bond.atom1 = atom1
        bond.atom2 = atom2
        bond.type = str(int(atom1.type) + n_bond_types)
        bond_num += 1
        bond.num = bond_num
        frame.insert_bond(bond)
        
frame.write_lmp_data(args.output[0])
