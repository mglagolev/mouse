#!/usr/bin/env python

import sys

def clusteringCore(aggregate, allAtomsCheckedInums, atom, start = False):
#Using atom numbers to for unique identification of atoms
	if start:
		aggregate.append(atom)
		allAtomsCheckedInums.append(atom.inum)
	neighbors = atom.neighbors
	neighborsToCheck = []
	for neighbor in neighbors:
		if neighbor.inum not in allAtomsCheckedInums:
			neighborsToCheck.append(neighbor)
			aggregate.append(neighbor)
			allAtomsCheckedInums.append(neighbor.inum)
	for nextatom in neighborsToCheck:
		clusteringCore(aggregate, allAtomsCheckedInums, nextatom)

def firstMissingInteger(somelist):		#find the first missing positive integer in a list of
	result = 1
	while 1:
		try:
			somelist.index(result)
			result += 1
		except ValueError: return result
		

def makeNeighborlistsFromBonds(config):
	for atom in config.atoms():
		atom.neighbors = []
	for bond in config.bonds():
		atom1 = bond.atom1
		atom2 = bond.atom2
		atom1.add_neighbor(atom2)
		atom2.add_neighbor(atom1)


def makeNeighborlistsFromDistances(config, threshold = 1.2):
	for atom in config.atoms():
		atom.neighbors = []
	atoms = config.atoms()
	for i in range(len(atoms)):
		atom1 = atoms[i]
		for atom2 in atoms[i+1:]:
			r = config.interatom_vector(atom1.num, atom2.num)[0].length()
			if r <= threshold:
				atom1.add_neighbor(atom2)
				atom2.add_neighbor(atom1)
			

def determineClusters(config):
	all_aggregates = []
	allAtomsCheckedInums = []
	startAtomInum = 1
	while startAtomInum <= config.n_atom():
		aggregate = []
		startatom = config.atom_by_inum(startAtomInum)
		clusteringCore(aggregate, allAtomsCheckedInums, startatom, start = True)
		aggregate.sort(key = lambda x: x.num)
		all_aggregates.append(aggregate)
		sys.stderr.write("\rChecked: " + str(len(allAtomsCheckedInums)))
		startAtomInum = firstMissingInteger(allAtomsCheckedInums)
	return all_aggregates


def assignMoleculesFromBonds(config):
	makeNeighborlistsFromBonds(config)
	molecules = determineClusters(config)
	for nmol in range(len(molecules)):
		for atom in molecules[nmol]:
			atom.mol_id = str(nmol)
