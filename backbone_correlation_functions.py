#!/usr/bin/env python

from lammpack_misc import *
from lammpack_types import *
import numpy
import sys

def CalculateBackboneCorrelations(config, bond_types, k_max = None):
	natom = config.n_atom()
	nbonds = config.n_bond()
	sys.stderr.write('Natom: '+str(natom)+' Nbonds: '+str(nbonds)+'\n')
	molecules_index, molecules_atoms = GetMoleculeAtomlists(config)
	all_bonds = []
	for i in range(len(molecules_index)):
		molecule_id = molecules_index[i]
		molecule_atoms = molecules_atoms[i]
		bonds = []
		for j in range(nbonds):
			bond = config.bond(j)
			try:
				bond_types.index(bond.type)
				try:
					molecule_atoms.index(bond.atom1.num)
					bonds.append(bond)
				except ValueError: pass
			except ValueError: pass
		bonds.sort(key = lambda bond: bond.atom1.num)
		all_bonds.append(bonds)
	molecules_nbonds = []
	for i in range(len(molecules_index)):
		molecule_nbonds = len(all_bonds[i])
		molecules_nbonds.append(molecule_nbonds)
	nbonds_min = min(molecules_nbonds)
	ck_sum = [0.] * nbonds_min
	ck_norm = [0] * nbonds_min
	for i in range(len(molecules_index)):
		bonds = all_bonds[i]
		for j in range(nbonds_min):
			if k_max is not None:
				k_end = min(j + 1 + int(k_max), nbonds_min)
			else:
				k_end = nbonds_min
			for k in range(j+1, k_end):
				bond1_vector, _ = config.bond_vector(bonds[j].num)
				bond2_vector, _ = config.bond_vector(bonds[k].num)
				bond1_vector_std = [bond1_vector.x, bond1_vector.y, bond1_vector.z]
				bond2_vector_std = [bond2_vector.x, bond2_vector.y, bond2_vector.z]
				ck_sum[k-j] += numpy.dot(bond1_vector_std, bond2_vector_std)/numpy.linalg.norm(bond1_vector_std)/numpy.linalg.norm(bond2_vector_std)
				ck_norm[k-j] += 1
	for i in range(len(ck_sum)):
		if ck_norm[i] > 0:
			print i, float(ck_sum[i]) / float(ck_norm[i])
		
