#!/usr/bin/env python3

from lammpack_misc import *
from lammpack_types import *
import vector3d
import numpy as np
import sys
import pdb

def hash8(s):
	return abs(hash(s)) % (10 ** 8)


def createNumpyBondsArrayFromConfig(config, allowedBondTypes = [], allowedResTypes = [], allowedMolIds = [], allowedBondNums = []):
	"""Create numpy 2d array, containing subarrays with bond numbers, bond vector coordinates, and bond center coordinates.
	   If the allowed bond types, res_types, mol_ids or numbers are specified, those should be in the corresponding lists.
	   The properties of the atoms should be in the lists for both of the atoms.
	   Returns: np.array((bond numbers, molecule ids, bond vector x, bond vector y, bond vector z,
	   bond center x, bond center y, bond center z))"""
	resTypesStr, molIdsStr = [], []
	bondNums, resTypes, molIds = [], [], []
	bx, by, bz = [], [], []
	rx, ry, rz = [], [], []
	for bond in config.bonds():
		if ( allowedBondTypes is None or bond.type in allowedBondTypes or len(allowedBondTypes) == 0 ) \
		and ( allowedResTypes is None or ( bond.atom1.res_type in allowedResTypes and bond.atom2.res_type in allowedResTypes ) or len(allowedResTypes) == 0 ) \
		and ( allowedMolIds is None or ( bond.atom1.mol_id in allowedMolIds and bond.atom2.mol_id in allowedMolIds ) or len(allowedMolIds) == 0 ) \
		and ( allowedBondNums is None or bond.num in allowedBondNums or len(allowedBondNums) == 0):
			b, r = config.bond_vector_by_num(bond.num)
			bondNums.append(bond.num)
			resTypesStr.append(bond.atom1.res_type)
			molIdsStr.append(bond.atom1.mol_id)
			bx.append(b.x)
			by.append(b.y)
			bz.append(b.z)
			rx.append(r.x)
			ry.append(r.y)
			rz.append(r.z)
	resTypes = list(map(hash8, resTypesStr))
	molIds = list(map(hash8, molIdsStr))
	return np.array((np.array(bondNums), np.array(resTypes), np.array(molIds), np.array(bx), np.array(by), np.array(bz), np.array(rx), np.array(ry), np.array(rz)))


def createNumpyAtomsArrayFromConfig(config, allowedAtomTypes = [], allowedResTypes = [], allowedMolIds = [], allowedAtomNums = []):
	"""Create numpy 2d array, containing subarrays with atom numbers and atom coordinates.
	   Returns: np.array((atom numbers, atom res_types, atom mol_ids, atom coords x, atom coords y, atom coords z))"""
	resTypesStr, molIdsStr = [], []
	atomNums, resTypes, molIds = [], [], []
	x, y, z = [], [], []
	for atom in config.atoms():
		if ( atom.type in allowedAtomTypes or len(allowedAtomTypes) == 0 ) \
		and ( atom.res_type in allowedResTypes or len(allowedResTypes) == 0 ) \
		and (atom.mol_id in allowedMolIds or len(allowedMolIds) == 0) \
		and (atom.num in allowedAtomNums or len(allowedAtomNums) == 0):
			atomNums.append(atom.num)
			resTypesStr.append(atom.res_type)
			molIdsStr.append(atom.mol_id)
			x.append(atom.pos.x)
			y.append(atom.pos.y)
			z.append(atom.pos.z)
	resTypes = map(hash8, resTypesStr)
	molIds = map(hash8, molIdsStr)
	return np.array((atomNums, resTypes, molIds, np.array(x), np.array(y), np.array(z)))


def calculateRdfForReference(config, refAtom, npAtoms, rmin = 0., rmax = -1., nbin = 100, excludeSelf = True, sameMolecule = True):
	"""Returns numpy histogram for distribution of distance between the atom positions in npAtoms array and the reference atom refAtom."""
	if rmax < 0.: rmax = min(config.box().x, config.box().y, config.box().z) / 2.
	drx = np.add(npAtoms[3], -1.* refAtom.pos.x)
	dry = np.add(npAtoms[4], -1.* refAtom.pos.y)
	drz = np.add(npAtoms[5], -1.*  refAtom.pos.z)
	drxTrim = np.add(drx, -1. * config.box().x * np.around(drx / config.box().x))
	dryTrim = np.add(dry, -1. * config.box().y * np.around(dry / config.box().y))
	drzTrim = np.add(drz, -1. * config.box().z * np.around(drz / config.box().z))
	if excludeSelf: notSelf = npAtoms[0] != refAtom.num
	else: notSelf = 1.
	if not sameMolecule:
		emptyMolId = npAtoms[2] == hash8('')
		if np.sum(emptyMolId) > 0:
			raise NameError("We should omit the atoms belonging to the same molecules, but some mol_ids are empty")
		notSameMolecule = npAtoms[2] != hash8(refAtom.mol_id)
	else: notSameMolecule = 1.
	drTrim = notSelf * notSameMolecule * np.sqrt(drxTrim**2 + dryTrim**2 + drzTrim**2)
	drTrimMasked = np.ma.masked_equal(drTrim, 0.)
	if drTrimMasked.count() > 0:
		return np.histogram(drTrimMasked.compressed(), bins = nbin, range = (rmin, rmax))
	else:
		raise NameError("No values to calculate")


def calculateAveCosSqForReference(config, refBond, npBonds, rcut, excludeSelf = True, sameMolecule = True):
	bRef, bRefCenter = config.bond_vector_by_num(refBond.num)
	if bRef.x != 0. or bRef.y != 0. or bRef.z != 0.:
		bx, by, bz = npBonds[3], npBonds[4], npBonds[5]
		drx = npBonds[6] - bRefCenter.x
		dry = npBonds[7] - bRefCenter.y
		drz = npBonds[8] - bRefCenter.z
		drxTrim = drx - config.box().x * np.around(drx / config.box().x)
		dryTrim = dry - config.box().y * np.around(dry / config.box().y)
		drzTrim = drz - config.box().z * np.around(drz / config.box().z)
		rsq = drxTrim**2 + dryTrim**2 + drzTrim**2
		if rcut >= 0.: inRange = rsq <= rcut**2
		else: inRange = 1.
		if excludeSelf: notSelf = npBonds[0] != refBond.num
		else: notSelf = 1.
		if not sameMolecule:
			emptyMolId = npBonds[2] == hash8('')
			if np.sum(emptyMolId) > 0:
				raise NameError("We should omit the atoms belonging to the same molecules, but some mol_ids are empty")
			notSameMolecule = npBonds[2] != hash8(refBond.atom1.mol_id)
		else: notSameMolecule = 1.
		cosSqNormed = inRange * notSelf * notSameMolecule * ( bRef.x * bx + bRef.y * by + bRef.z * bz )**2 / ( bRef.x**2 + bRef.y**2 + bRef.z**2 ) / ( bx**2 + by**2 + bz**2)
		cosSqNormedMasked = np.ma.masked_equal(cosSqNormed, 0.)
		if cosSqNormedMasked.count() > 0:
			return np.ma.average(cosSqNormedMasked)
		else:
			raise NameError("No values to calculate")
	else:
		raise NameError("Zero length of reference bond " + str(refBond.num))


def CalculatePairCorrelationFunctions(config, atomtypes = [], rmin = 0., rmax = 0., nbin = 100, sameMolecule = True):
	if rmax < 0.:
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
	types.append('AAAAA')
	atom_counter['AAAAA'] = 0
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
			rdfdata[pair_string] = np.ndarray((nbin))
			norm[pair_string] = 0
	rdfs['data'] = rdfdata
	rdfs['norm'] = norm
	npAtomsByType = {}
	for atype in types:
		npAtomsByType[atype] = createNumpyAtomsArrayFromConfig(config, allowedAtomTypes = [atype])
	for atom in atoms:
		type1 = atom.type
		if type1 in types:
			for type2 in types:
				pair_types = [type1, type2]
				pair_types.sort()
				pair_string = str(pair_types[0]) + '_' + str(pair_types[1])
				try:
					rdf = calculateRdfForReference(config, atom, npAtomsByType[type2], rmin = rmin, rmax = rmax, nbin = nbin, sameMolecule = sameMolecule)
					rdfs['data'][pair_string] = np.add(rdfs['data'][pair_string],rdf[0])
					rdfs['norm'][pair_string] += np.sum(rdf[0])
				except NameError: pass
	for pair_string in rdfs['norm']:
		rdfs['norm'][pair_string] = float(rdfs['norm'][pair_string]) / 2. / config.box().x / config.box().y / config.box().z
	sys.stderr.write(str(rdfs))
	return rdfs


def CalculateOrientationOrderParameter(config, bondtypes, rmin = 0., rmax = 0., storeAsAtomtypes = False, mode = 'average', smin = -0.5, smax = 1.0, sbins = 75, referenceResTypes = [], sameMolecule = True):
	nbonds = config.n_bond()
	cos_sq = 0.
	i_s = 0
	progress = 0
	if rmax == 0.:
		rmax = config.box().length()/2.
	if mode == 'histo':
		s_hist = createHistogram(smin, smax, sbins)
		s_hist_norm_total = 0.
		for i in range(sbins):
			sleft, sright = smin + i * s_hist["step"], smin + (i + 1) * s_hist["step"]
			costhetaleft, costhetaright = (2./3. * (sleft + 0.5))**0.5, (2./3. * (sright + 0.5))**0.5
			s_hist["norm"][i] = 1. / (costhetaleft + costhetaright)
	if storeAsAtomtypes:
		atomOrderingList = []
	npBonds = createNumpyBondsArrayFromConfig(config, allowedBondTypes = bondtypes)
	for bond in config.bonds():
		if ( bondtypes is None or bond.type in bondtypes or len(bondtypes) == 0 ) \
		and ( referenceResTypes is None or ( bond.atom1.res_type in referenceResTypes and bond.atom2.res_type in referenceResTypes ) or len(referenceResTypes) == 0 ):
			try:
				cosSq = calculateAveCosSqForReference(config, bond, npBonds, rmax, excludeSelf = True, sameMolecule = sameMolecule)
				if mode == 'histo':
					updateHistogram(1.5 * cosSq - 0.5, s_hist)
				elif mode == 'average':
					cos_sq += cosSq
					i_s += 1
				if storeAsAtomtypes:
					bondAtomNums = [bond.atom1.num, bond.atom2.num]
					for atomNum in bondAtomNums:
						try:
							atomOrderingList[[x["num"] for x in atomOrderingList].index(atomNum)]["values"].append(1.5 * cosSq - 0.5)
						except ValueError:
							atomOrderingList.append({ "num" : atomNum, "values" : [1.5 * cosSq - 0.5] })
			except NameError: pass
	if storeAsAtomtypes:
		for atomOrdering in atomOrderingList:
			atom = config.atom_by_num(atomOrdering["num"])
			atomOrderingMean = np.mean(atomOrdering["values"])
			if atomOrderingMean >= 0.: atom.type = "C" + "%02d" % int(atomOrderingMean * 100)
			else: atom.type = "Cm" + "%02d" % int(atomOrderingMean * -100)
	if mode == 'average':
		if i_s > 0:
			s = 1.5 * cos_sq / i_s - 0.5
		return s
	elif mode == 'histo':
		normHistogram(s_hist)
		return s_hist


def CalculateCrystallinityParameter(config, bondtypes, reference_vector=vector3d.Vector3d(0. ,0. ,0. ), storeAsAtomtypes = False, mode = 'average', smin = -0.5, smax = 1.0, sbins = 75 ):
	if mode == 'histo':
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
		if (bondtypes is None) or (bond.type in bondtypes): bondlist.append(bond.num)
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
	if mode == 'average':
		npBonds = createNumpyBondsArrayFromConfig(config, allowedBondTypes = bondtypes)
		refBond = Bond()
		refBond.atom1, refBond.atom2 = Atom(), Atom()
		refBond.atom1.pos = vector3d.Vector3d(0., 0., 0.)
		refBond.atom2.pos = ave_b
		cosSq = calculateAveCosSqForReference(config, refBond, npBonds, -1.)
		return 1.5 * cosSq - 0.5
	elif mode == 'histo':
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
			updateHistogram(s_value, s_hist)
		if storeAsAtomtypes:
			for atomCryst in atomOrderingList:
				atom = config.atom_by_num(atomCryst["num"])
				atomCrystallinity = np.mean(atomCryst["values"])
				if atomCrystallinity >= 0.: atom.type = "C" + "%02d" % int(atomCrystallinity * 100)
				else: atom.type = "Cm" + "%02d" % int(atomCrystallinity * -100)
		normHistogram(s_hist, normInternalNorm = True)
		return s_hist
		
def CalculateComRdfs(config, rmin, rmax, nbin):
	""" Calculate radial distribution function for different atomtypes depending on the distance
	from the center of mass of the system """
	import ordering_functions as of
	import copy
	cpconf = copy.deepcopy(config)
	cpconf.convertCoordsCutToUncut()
	center = Atom()
	center.pos = cpconf.com()
	# Create a set of atom types
	atypes = set(list(atom.type for atom in cpconf.atoms())) # Create a set of unique atom types
	#Prepare rdf 
	rdfs = {}
	distances, volumes = [], []
	for i in range(nbin):
		distances.append(rmin + (rmax - rmin) * (i+0.5) / nbin)
		r1 = rmin + (rmax - rmin) * i / nbin
		r2 = rmin + (rmax - rmin) * (i + 1) / nbin
		volumes.append(4./3.*np.pi*(r2**3 - r1**3))
	rdfs['r'], rdfs['v'] = distances, volumes
	rdfdata = {}
	norm = {}
	for atype in atypes:
		npAtomsByType = of.createNumpyAtomsArrayFromConfig(cpconf, allowedAtomTypes = [atype])
		rdf = of.calculateRdfForReference(cpconf, center, npAtomsByType, rmin = rmin, rmax = rmax, nbin = nbin, excludeSelf = False, sameMolecule = True)
		rdfdata[atype] = rdf[0]
		norm[atype] = np.sum(rdf[0])
	rdfs['data'] = rdfdata
	rdfs['norm'] = norm
	return rdfs
