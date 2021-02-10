#!/usr/bin/env python3

import sys, argparse
from lammpack_types import *

parser = argparse.ArgumentParser(description = 'Replicate the .pdb file along the axes')

parser.add_argument('in_pdb', metavar = 'PDB_IN', type = str, nargs = '?', default = 'in.pdb', help = 'input .pdb file')

parser.add_argument('out_pdb', metavar = 'PDB_OUT', type = str, nargs = '?', default = 'out.pdb', help = 'output .pdb file')

parser.add_argument('--x', type = int, nargs = '?', default = 1, help = 'number of copies along the X axis')

parser.add_argument('--y', type = int, nargs = '?', default = 1, help = 'number of copies along the Y axis')

parser.add_argument('--z', type = int, nargs = '?', default = 1, help = 'number of copies along the Z axis')

args = parser.parse_args()

replicate_x = range(int(-1 * (args.x / 2)),int(args.x - args.x / 2))
replicate_y = range(int(-1 * (args.y / 2)), int(args.y - args.y / 2))
replicate_z = range(int(-1 * (args.z / 2)), int(args.z - args.z / 2))

frame = Config()

frame.read_pdb(args.in_pdb)
box = frame.box()
atoms = frame.atoms()[:]
natom = len(atoms)
bonds = frame.bonds()[:]
i_replica = 0
for i in replicate_x:
	for j in replicate_y:
		for k in replicate_z:
			if i!=0 or j!=0 or k!=0:
				i_replica += 1
				for atom in atoms:  
					newatom = Atom()
					newatom.pos.x = atom.pos.x + i * box.x
					newatom.pos.y = atom.pos.y + j * box.y
					newatom.pos.z = atom.pos.z + k * box.z
					newatom.num = atom.num + i_replica * natom
					newatom.type = atom.type
					newatom.element = atom.element
					newatom.id = atom.id
					newatom.mol_id = atom.mol_id
					newatom.res_type = atom.res_type
					frame.insert_atom(newatom)
				for bond in bonds:
					newbond = Bond()
					newbond.atom1.num = bond.atom1.num + i_replica * natom
					newbond.atom2.num = bond.atom2.num + i_replica * natom
					frame.insert_bond(newbond)
frame.write_pdb(args.out_pdb, False)
