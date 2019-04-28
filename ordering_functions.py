#!/usr/bin/env python

from lammpack_misc import *
from lammpack_types import *
import vector3d
import numpy as np
import sys
import pdb

def CalculatePairCorrelationFunctions(config, atomtypes = [], rmin = 0., rmax = 0., nbin = 100, sameMolecule = True):
	if rmax <= 0.:
		rmax = min(config.box().x, config.box().y, config.box().z) / 2.
	types = []
	atoms = config.atoms()
	atom_counter = {}
	for atom in atoms:
		if atomtypes is None or len(atomtypes) == 0 or ( atom.type in atomtypes ):
			try:
				types.index(atom.type)
				atom_counter[atom.type] += 1
			except ValueError:
				types.append(atom.type)
				atom_counter[atom.type] = 1
	types.sort()
	rdfs = {}
	distances = []
	for i in range(nbin):
		distances.append(rmin + (rmax - rmin) * (i+0.5) / nbin)
	rdfs['r'] = distances	
	volumes = []
	for i in range(nbin):
		r1 = rmin + (rmax - rmin) * i / nbin
		r2 = rmin + (rmax - rmin) * (i + 1) / nbin
		volumes.append(4./3.*np.pi*(r2**3 - r1**3))
	rdfs['v'] = volumes
	rdfdata = {}
	norm = {}
	for i in range(len(types)):
		for j in range(i,len(types)):
			pair_string = str(types[i]) + '_' + str(types[j])
			rdfdata[pair_string] = [0] * nbin
			norm[pair_string] = 0
	rdfs['norm'] = norm
	for i in range(len(atoms)):
		for j in range(i+1, len(atoms)):
			type1 = atoms[i].type
			type2 = atoms[j].type
			if (type1 in types) and (type2 in types):
				mol1 = atoms[i].mol_id
				mol2 = atoms[j].mol_id
				if sameMolecule or (mol1 != mol2):
					pair_types = [type1, type2]
					pair_types.sort()
					pair_string = str(pair_types[0]) + '_' + str(pair_types[1])
					rdfs['norm'][pair_string] += 1
					dist = vector_pbc_trim(atoms[i].pos - atoms[j].pos, config.box()).length()
					ihist = int( nbin * (dist - rmin) / (rmax - rmin))
					if ihist >= 0 and ihist < nbin:
						rdfdata[pair_string][ihist] += 1
	rdfs['data'] = rdfdata
	for pair_string in rdfs['norm']:
		rdfs['norm'][pair_string] = float(rdfs['norm'][pair_string]) / config.box().x / config.box().y / config.box().z
	return rdfs

def CalculateOrientationOrderParameter(config, bondtypes, rmin = 0., rmax = 0., mode = 'average', smin = -0.5, smax = 1.5, sbins = 75):
	nbonds = config.n_bond()
	cos_sq = 0.	
	i_s = 0
	progress = 0
	bondlist = []
	if rmax == 0.:
		rmax = config.box().length()/2.
	if mode == 'histo':
		s_hist = createHistogram(smin, smax, sbins)
	for bond in config.bonds():
		try:
			bondtypes.index(bond.type)
			bondlist.append(bond.num)
		except: pass
	for i in range(len(bondlist)):
		for j in range(i+1, len(bondlist)):
			b1, r1 = config.bond_vector_by_num(bondlist[i])
			b2, r2 = config.bond_vector_by_num(bondlist[j])
			d = vector3d.Vector3d()
			d = r2.__sub__(r1)
			d.x = d.x - config.box().x * int(round(d.x / config.box().x))
			if abs(d.x) <= rmax:
				d.y = d.y - config.box().y * int(round(d.y / config.box().y))
				if abs(d.y) <= rmax:
					d.z = d.z - config.box().z * int(round(d.z / config.box().z))
					if abs(d.z) <= rmax:
						if d.length() <= rmax and d.length > rmin:
							cos_sq_value = (b1.x*b2.x+b1.y*b2.y+b1.z*b2.z)**2/b1.length_sq()/b2.length_sq()
							if mode == 'average':
								cos_sq += cos_sq_value
								i_s += 1
							elif mode == 'histo':
								updateHist(1.5 * cos_sq_value - 0.5, s_hist)
			progress += 1
			if progress % 1000 == 0:
				sys.stderr.write("\r"+str(progress)+" "+str((1.5 * cos_sq)/float(i_s)-0.5)+" "+str(i_s))
	if mode == 'average':
		if i_s > 0:
			s = 1.5 * cos_sq / i_s - 0.5
		return s
	elif mode == 'histo':
		normHistogram(s_hist)
	 	return s_hist

def CalculateCrystallinityParameter(config, bondtypes, reference_vector=vector3d.Vector3d(0. ,0. ,0. ), storeAsAtomtypes = False, mode = 'average', smin = -0.5, smax = 1.0, sbins = 75 ):
	if mode == 'average':
		s = 0.
		i_s = 0
	elif mode == 'histo':
		s_hist = createHistogram(smin, smax, sbins)
		s_hist_norm_total = 0.
		for i in range(sbins):
			sleft, sright = smin + i * s_hist["step"], smin + (i + 1) * s_hist["step"]
			costhetaleft, costhetaright = (2./3. * (sleft + 0.5))**0.5, (2./3. * (sright + 0.5))**0.5
			s_hist["norm"][i] = 1. / (costhetaleft + costhetaright)
	progress = 0
	bondlist = []
	if storeAsAtomtypes:
		atomOrderingList = []
	for bond in config.bonds():
		if bondtypes is None: bondlist.append(bond.num)
		try:
			bondtypes.index(bond.type)
			bondlist.append(bond.num)
		except: pass
	if reference_vector.length_sq() == 0.:
		ave_b = vector3d.Vector3d()
		i_b = 0
		for i in range(len(bondlist)):
			b1, r1 = config.bond_vector_by_num(bondlist[i])
			ave_b = ave_b + b1
			i_b += 1
		if i_b > 0:
			ave_b.scale(1./float(i_b))
	else:
		ave_b = reference_vector
	for i in range(len(bondlist)):
		b1, r1 = config.bond_vector_by_num(bondlist[i])
		s_value = 1.5 * (b1.x*ave_b.x+b1.y*ave_b.y+b1.z*ave_b.z)**2/b1.length_sq()/ave_b.length_sq() - 0.5
		if storeAsAtomtypes:
			bondAtomNums = [config.bond_by_num(bondlist[i]).atom1.num, config.bond_by_num(bondlist[i]).atom2.num]
			for atomNum in bondAtomNums:
				try:
					atomOrderingList[[x["num"] for x in atomOrderingList].index(atomNum)]["values"].append(s_value)
				except ValueError:
					atomOrderingList.append({ "num" : atomNum, "values" : [s_value] })
		if mode == 'average':
			s += s_value
			i_s += 1
		elif mode == 'histo':
			updateHistogram(s_value, s_hist)
	if storeAsAtomtypes:
		for atomCryst in atomOrderingList:
			atom = config.atom_by_num(atomCryst["num"])
			atomCrystallinity = np.mean(atomCryst["values"])
			if atomCrystallinity >= 0.: atom.type = "C" + "%02d" % int(atomCrystallinity * 100)
			else: atom.type = "Cm" + "%02d" % int(atomCrystallinity * -100)
	if mode == 'average':
		if i_s > 0:
			s /= i_s
			return s
	elif mode == 'histo':
		normHistogram(s_hist, normInternalNorm = True)
		return s_hist
