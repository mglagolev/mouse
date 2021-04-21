#!/usr/bin/python3

from lammpack_misc import *
from lammpack_types import *
import numpy as np
import sys
import math
import vector3d
import copy
import operator
import re

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def Cart2Sph(vector):
	x = vector[0]
	y = vector[1]
	z = vector[2]
	rxy = np.sqrt(x**2 + y**2)
	r = np.sqrt(x**2 + y**2 + z**2)
	phi = np.arctan2(y, x)
	theta = np.arctan2(z, rxy)
	return [r, theta, phi]

def Vector3d2NpArray(vector):
	return np.array([vector.x, vector.y, vector.z])

def NpArray2Vector3d(array):
	return vector3d.Vector3d(array[0], array[1], array[2])

#Calculate the coordinates relative to position of the "key0" atom
def place_at_0(key0, coords):
	coords_key0 = coords[key0]
	for key in coords:
		coords[key] = coords[key] - coords_key0

#Rotate around the vector which goes through O1 in the XY plane perpendicular to the O1-A vector so that A is in the XY plane
def rotate_into_XY(key0, target, coords):
	rotation_vector = np.cross(coords[target] - coords[key0], np.array([0.,0.,1.]))
	rotation_angle = -1. * np.arcsin(coords[target][2] / np.linalg.norm( coords[target] - coords[key0]))
	rotate_group_around(key0, coords, rotation_vector, rotation_angle)

#Rotate around the vector which goes through O1 along Z axis so that O1-A is parallel to X axis
def rotate_into_X(key0, target, coords):
	rotation_vector = np.array([0., 0., 1.])
	rotation_angle = -1. * np.arctan2(coords[target][1], coords[target][0])
	rotate_group_around(key0, coords, rotation_vector, rotation_angle)

def rotate_group_around(rotation_key, coords, rotation_vector, rotation_angle):
	coords_rotation_key = coords[rotation_key]
	for key in coords:
		coords[key] = coords_rotation_key + (np.dot(rotation_matrix(rotation_vector,rotation_angle), coords[key] - coords_rotation_key))

def scale_group(coords, factor):
	for key in coords:
		coords[key] = coords[key] * factor

def scale_uniform(key0, coords, factor):
	coords_key0 = coords[key0]
	for key in coords:
		a = coords[key] - coords_key0
		coords[key] = coords_key0 + factor * a

def scale_along_vector(key0, scaling_vector, coords, factor):
	coords_key0 = coords[key0]
	sc_vec_norm = scaling_vector / np.linalg.norm(scaling_vector)
	for key in coords:
		a = coords[key] - coords_key0
		a_par = sc_vec_norm * np.dot(a, sc_vec_norm)
		a_ort = a - a_par
		coords[key] = coords_key0 + factor * a_par + a_ort

def translate_group(coords, translation_vector):
	for key in coords:
		coords[key] = coords[key] + translation_vector

def PrepareCoords(unit_type, fname, key):
	#This function should read raw coordinates from the relevant .dat file and make initial alignment (O1-A along X axis, B in XY plane)
	s = open(fname, 'r').read()
	coords = eval(s)
	return RotateCoords(unit_type, key, coords)

def RotateCoords(unitType, key, coords):
	if unitType == "SLS": key0, target = "A", "A_next"
	elif unitType == "B": key0, target = key, "C6"
	else: key0, target = key, "A"
	place_at_0(key0, coords)
	rotate_into_XY(key0, target, coords)
	rotate_into_X(key0, target, coords)
	return coords

def calculateUnitLength(unitType, key, coords):
	if unitType == "SLS": key0, target = "A", "A_next"
	else: key0, target = key, "A"
	return np.linalg.norm(coords[target] - coords[key0])

def FetchCoords(templatesList, key):
	templatesDict = {}
	for filename in templatesList:
		s = open(filename, 'r').read()
		lines = s.split('\n')
		l = lines[0].split()
		if l[0] == '#' and l[1] == "unit_type":
			unitType = l[2]
		coords = RotateCoords(unitType, key, eval(s))
		unitLength = calculateUnitLength(unitType, key, coords)
		if unitType in templatesDict:
			templatesDict[unitType].append([unitLength, coords])
		else:
			templatesDict[unitType] = []
			templatesDict[unitType].append([unitLength, coords])
	for unitType in templatesDict:
		templatesDict[unitType].sort(key = operator.itemgetter(0))
	for key in templatesDict:
		print key
	return templatesDict

def AlignB(configBox, coords, A_pos, B_pos, rotation_center, rotation_vector):
	R_norm = rotation_vector / np.linalg.norm(rotation_vector)
	AB_templ = Vector3d2NpArray(vector_pbc_trim(NpArray2Vector3d(coords["B"]) - NpArray2Vector3d(coords["A"]), configBox))
	AB_templ_norm = AB_templ / np.linalg.norm(AB_templ)
	AB_config = Vector3d2NpArray(vector_pbc_trim(B_pos - A_pos, configBox))
	AB_config_norm = AB_config / np.linalg.norm(AB_config)
	AB_templ_norm_x_R = np.cross(AB_templ_norm, R_norm)
	AB_config_norm_x_R = np.cross(AB_config_norm, R_norm)
	rotation_sin = np.dot(np.cross(AB_config_norm_x_R, AB_templ_norm_x_R), R_norm)
	rotation_cos = np.dot(AB_config_norm_x_R, AB_templ_norm_x_R)
	angle = -1 * np.arctan2(rotation_sin, rotation_cos)
	rotate_group_around(rotation_center, coords, rotation_vector, angle)
	return angle
	
def getClosestCoords(coordsList, CG_dist):
	for i in range(len(coordsList)):
		AA_dist = coordsList[i][0]
		if AA_dist > CG_dist:
			return coordsList[i][1]
	return coordsList[-1][1]

def Restore(scalingWeight, restore_type, configBox, template, cgAnchor_pos, cgTarget_pos, pendant_pos = None):
	"""
		The function restores the atomistic structure of the fragment using the template.
		The alignment of the template is based on the position of the anchor particle ('cgAnchor_pos'), the position of the target particle ('cgTarget_pos') and position of the side pendant.
		The particles of the template which should be superimposed upon anchor and target atom of the config are encoded in the 'restore_type' parameter (restore_type = ANCHOR-TARGET).
		The transformations include:
			- rotating the template so that the orientation of vector from anchor atom to target atom in the template and in the system coincide
			- scaling the template along the abovementioned vector with a weight: 0% - no scaling; 100% - the target atoms in the template and in the system fully coincide
			- placing the anchor atom of the template in the position of anchor atom of the existing configuration
			- rotation of the template around the same vector to determine the optimal position for the side pendant
	"""
	key0, target = restore_type.split("-")[0], restore_type.split("-")[1]
	sph = Cart2Sph(Vector3d2NpArray(vector_pbc_trim(cgTarget_pos - cgAnchor_pos, configBox)))
	CG_dist = sph[0]
	_coords = template.copy()
	rotate_group_around(key0, _coords, np.array([0.,1.,0.]), -1. * sph[1])
	rotate_group_around(key0, _coords, np.array([0.,0.,1.]), 1. * sph[2])
	AA_dist = np.linalg.norm(_coords[target] - _coords[key0])
	factor = CG_dist / AA_dist
	if factor <= 1.: weightedFactor = 1. + (factor - 1.) * scalingWeight[0]
	else:  weightedFactor = 1. + (factor - 1.) * scalingWeight[1]
	scale_along_vector(key0, _coords[target] - _coords[key0], _coords, weightedFactor)
	#scale_uniform(key0, _coords, weightedFactor)
	translate_group(_coords, Vector3d2NpArray(cgAnchor_pos) - _coords[key0])
	angle = 0.
	if pendant_pos is not None:
		rotation_vector = _coords[target] - _coords[key0]
		angle = AlignB(configBox, _coords, cgTarget_pos, pendant_pos, key0, rotation_vector)
	return _coords, weightedFactor, angle

def MakeAtomsFromCoords(coords, res_num, res_type, natom, cg_keys, substitute_keys):
	#Make a dictionary of atom objects from indexed coordinates
	atoms = []
	i = 1
	sequence = ["O1", "O1_SLS", "O1_SLE", "C2", "C3", "O4", "H5", "H10", "O10", "H11", "C6", "H7", "H8", "H9", "O1_next", "C3_prev", "C3_next", "A", "B", "A_next"]
	for key in sequence:
		try:
			position = coords[key]
			try:
				key = substitute_keys[key]
			except KeyError: pass
			atom = Atom()
			atom.res_num, atom.res_type, atom.type, atom.num = res_num, res_type, key, natom + i
			try: cg_keys.index(key)
			except ValueError: i += 1
			atom.pos = NpArray2Vector3d(position)
			atoms.append(atom)
		except KeyError: pass
	return atoms

def InsertAtoms(config, atoms, cg_keys):
	for atom in atoms:
		atomname = atom.type
		try:
			cg_keys.index(atomname)
			continue
		except ValueError:
			config.insert_atom(atom)

def getKeyAtom(atoms, key):
	for atom in atoms:
		if atom.type == key:
			return atom

def adjustPendant(aa_coords, configBox, cgPendant, sideScalingWeight, sideScalingData):
	aa_side_coords = { "C3" : aa_coords["C3"], "C6" : aa_coords["C6"], "H7" : aa_coords["H7"], "H8" : aa_coords["H8"], "H9" : aa_coords["H9"] }
	aa_side_coords = RotateCoords("B", "C3", aa_side_coords)
	aa_corrected_side_coords, sideScalingFactor, _ = Restore(sideScalingWeight, "C3-C6", configBox, aa_side_coords, NpArray2Vector3d(aa_coords["C3"]), cgPendant.pos)
	for sideKey in aa_corrected_side_coords: aa_coords[sideKey] = aa_corrected_side_coords[sideKey]
	UpdateScalingData(sideScalingData, cgPendant, sideScalingFactor)

def UpdateScalingData(old_data, particle, scalingFactor):
	try:
		if scalingFactor < old_data[0]:
			old_data[0], old_data[1] = scalingFactor, [particle.mol_id, particle.num]
		elif scalingFactor > old_data[2]:
			old_data[2], old_data[3] = scalingFactor, [particle.mol_id, particle.num]
	except IndexError:
		old_data.append(scalingFactor)
		old_data.append([particle.mol_id, particle.num])
		old_data.append(scalingFactor)
		old_data.append([particle.mol_id, particle.num])

#def pla_reverse_map(config, sls_fname, slm_fname, sle_fname):
def pla_reverse_map(config, templatesList, scalingWeight, sideScalingWeight, doAdjustPendant):
	key, key_next = "C3_prev", "C3"
	cg_keys = [ "A", "A_next", "B", "C3_prev", "O1_next" ]
	substitute_keys = {}
	#key, key_next = "O1", "O1"
	#cg_keys = [ "A", "A_next", "B", "C3_prev", "O1_this" ]
	#substitute_keys = { "O1" : "O1_this", "O1_next" : "O1", "O1_SLS" : "O1" }
	templatesDict = FetchCoords(templatesList, key)
	aa_config = Config()
	scalingData, sideScalingData  = [], []
	aa_config.title = config.title
	aa_config.set_box(config.box())
	molecules_index, molecules_atoms = GetMoleculeAtomlists(config)
	for i in range(len(molecules_index)):
		mol_id = molecules_index[i]
		atoms_list = molecules_atoms[i]
		atoms_list.sort()
		npoly = len(atoms_list) / 2
		#Position SLS residue based on first two A beads
		cg_A1 = config.atom_by_num(atoms_list[0])
		cg_B1 = config.atom_by_num(atoms_list[1])
		cg_A2 = config.atom_by_num(atoms_list[2])
		if cg_A1.type != "A1" or cg_A2.type != "A2" or cg_B1.type != "B1":
			raise NameError("The types of atoms #1 #2 and #3 in molecule "+str(mol_id)+" are not A1, B1, A2")
		CG_dist = vector_pbc_trim(cg_A2.pos - cg_A1.pos, config.box()).length()
		template = getClosestCoords(templatesDict["SLS"], CG_dist)
		aa_coords, scalingFactor, rotationAngle = Restore(scalingWeight, "A-A_next", config.box(), template, cg_A1.pos, cg_A2.pos, cg_B1.pos)
		print "Unit ", cg_A1.mol_id, cg_A1.num, "rotation:", rotationAngle
		UpdateScalingData(scalingData, cg_A2, scalingFactor)
		if doAdjustPendant == True: adjustPendant(aa_coords, config.box(), cg_B1, sideScalingWeight, sideScalingData)
		aa_atoms = MakeAtomsFromCoords(aa_coords, mol_id, 'CG'+str(npoly)+'LA', aa_config.n_atom(), cg_keys, substitute_keys)
		InsertAtoms(aa_config, aa_atoms, cg_keys)
		aa_key = getKeyAtom(aa_atoms, key_next)
		for j in range(1,npoly-1):
			cg_A = config.atom_by_num(atoms_list[j*2])
			cg_B = config.atom_by_num(atoms_list[j*2+1])
			if cg_A.type[0] != "A" or cg_B.type[0] != "B":
				raise NameError("The types of atoms #1, #2 in unit " + str(j+1) + " of molecule " + str(mol_id) + " are not A*, B*")
			CG_dist = vector_pbc_trim(cg_A.pos - aa_key.pos, config.box()).length()
			template = getClosestCoords(templatesDict["SLM"], CG_dist)
			aa_coords, scalingFactor, rotationAngle = Restore(scalingWeight, key+"-A", config.box(), template, aa_key.pos, cg_A.pos, cg_B.pos)
			print "Unit ", cg_A1.mol_id, cg_A1.num, "rotation:", rotationAngle
			UpdateScalingData(scalingData, cg_A, scalingFactor)
			if doAdjustPendant: adjustPendant(aa_coords, config.box(), cg_B, sideScalingWeight, sideScalingData)
			aa_atoms = MakeAtomsFromCoords(aa_coords, mol_id, 'CG'+str(npoly)+'LA', aa_config.n_atom(), cg_keys, substitute_keys)
			InsertAtoms(aa_config, aa_atoms, cg_keys)
			aa_key = getKeyAtom(aa_atoms, key_next)
		cg_A = config.atom_by_num(atoms_list[-2])
		cg_B = config.atom_by_num(atoms_list[-1])
		if cg_A.type[0] != "A" or cg_B.type[0] != "B":
			raise NameError("The types of atoms #1, #2 in unit " + str(npoly) + " of molecule " + str(mol_id) + " are not A*, B*")
		CG_dist = vector_pbc_trim(cg_A.pos - aa_key.pos, config.box()).length()
		template = getClosestCoords(templatesDict["SLE"], CG_dist)
		aa_coords, scalingFactor, rotationAngle = Restore(scalingWeight, key+"-A", config.box(), template, aa_key.pos, cg_A.pos, cg_B.pos)
		print "Unit ", cg_A1.mol_id, cg_A1.num, "rotation:", rotationAngle
		UpdateScalingData(scalingData, cg_A, scalingFactor)
		if doAdjustPendant == True: adjustPendant(aa_coords, config.box(), cg_B, sideScalingWeight, sideScalingData)
		aa_atoms = MakeAtomsFromCoords(aa_coords, mol_id, 'CG'+str(npoly)+'LA', aa_config.n_atom(), cg_keys, substitute_keys)
		InsertAtoms(aa_config, aa_atoms, cg_keys)
		scaling_output = "\rProcessing molecule "+str(mol_id)+" Scaling min.: {:.3f}".format(scalingData[0])+" max.: {:.3f}".format(scalingData[2])
		if doAdjustPendant:
			scaling_output += " side scaling min.: {:.3f}".format(sideScalingData[0])+" max.: {:.3f}".format(sideScalingData[2])
		sys.stderr.write(scaling_output)
	sys.stdout.write("\nMinimum scaling factor: " + str(scalingData[0]) + " (molecule " + str(scalingData[1][0]) + ", atom " + str(scalingData[1][1]) + ")\n")
	sys.stdout.write("Maximum scaling factor: " + str(scalingData[2]) + " (molecule " + str(scalingData[3][0]) + ", atom " + str(scalingData[3][1]) + ")\n")
	if doAdjustPendant:
		sys.stdout.write("\nMinimum side scaling factor: " + str(sideScalingData[0]) + " (molecule " + str(sideScalingData[1][0]) + ", atom " + str(sideScalingData[1][1]) + ")\n")
		sys.stdout.write("Maximum side scaling factor: " + str(sideScalingData[2]) + " (molecule " + str(sideScalingData[3][0]) + ", atom " + str(sideScalingData[3][1]) + ")\n")
	return aa_config

def findAtom(config, atomType, atomResNum, startIndex = 0, direction = "forward"):
	if direction == "forward":
		atoms = config.atoms()[startIndex:]
	elif direction == "backward":
		atoms = list(reversed(config.atoms()[:startIndex]))
	for atom in atoms:
		if atom.type == atomType and atom.res_num == atomResNum:
			return atom

def writeTemplate(template, unitType, unitLength, fileName):
	out = open(fileName, 'wb')
	out.write('# ' + 'unit_type ' +  unitType + '\n')
	out.write('# ' + 'length ' + str(unitLength) + '\n')
	out.write('{\n')
	for key in template:
		out.write("\"" + key + "\" : np.array([ " + str(template[key].pos.x) + " , " + str(template[key].pos.y) + " , " + str(template[key].pos.z) + " ]),\n")
	out.write('}')
	out.close()

def writeUnitStructures(config, unitType, fileName):
	natoms = len(config.atoms())
	atomNumOffset = 0
	templateCounter = 0	
	while atomNumOffset < natoms:
		try:
			template, unitVector, maxAtom = getUnitStructure(config, unitType, atomNumOffset)
			atomNumOffset = maxAtom + 1
			templateCounter += 1
			writeTemplate(template, unitType, unitVector.length(), fileName + str(templateCounter))
		except: break

def getUnitStructure(config, unitType, atomNumOffset = 0):
#Determine variables:
	keyAtom = "C3"
#Determine set of required atom types
	lists = { "SLS" : ["O1", "C2", "C3", "O4", "H5", "H10", "C6", "H7", "H8", "H9"], \
	"SLM" : ["O1", "C2", "C3", "O4", "H5", "C6", "H7", "H8", "H9"], \
	"SLE" : ["O1", "C2", "C3", "O4", "H5", "O10", "H11", "C6", "H7", "H8", "H9"] }
	A_lists = { "SLS" : {"O1_next":15.9994, "O1":15.9994, "C2":12.011, "C3":12.011, "O4":15.9994, "H5":1.008, "H10":1.008}, \
	"SLM" : {"O1_next":15.9994, "C2":12.011, "C3":12.011, "O4":15.9994, "H5":1.008}, \
	"SLE" : {"C2":12.011, "C3":12.011, "O4":15.9994, "H5":1.008, "O10":15.9994, "H11":1.008} }
	B_lists = { "SLS" : {"C6":12.011, "H7":1.008, "H8":1.008, "H9":1.008}, \
	"SLM" : {"C6":12.011, "H7":1.008, "H8":1.008, "H9":1.008}, \
	"SLE" : {"C6":12.011, "H7":1.008, "H8":1.008, "H9":1.008} }
#Sort configuration atoms by number
	config.trim()
#Available atom identifiers:
#Monomer unit number (1, 2, 3), "res_num"; monomer unit type (SLS, SLM, SLE) "res_type"; atom type (O1...H9), "type"; atom number (enumeration across all molecules), "num"
	template_counter = 0
	maxAtomNum = 0
	for i in range(len(config.atoms())):
		atom = config.atoms()[i]
		if atom.num >= atomNumOffset:
			resNum = atom.res_num
			if atom.res_type == unitType and atom.type == lists[unitType][0]:
				template = {}
				for requiredType in lists[unitType]:	#Adding all of the atoms of the atomistic structure
					template[requiredType] = findAtom(config, requiredType, resNum, i, "forward")
				if unitType != "SLE":	#Adding O1 atom of the next monomer unit (to determine position of A bead)
					template["O1_next"] = findAtom(config, "O1", resNum + 1, i, "forward")
				if unitType != "SLS": #Adding C3 atom of the previous monomer unit (required for the template)
					template["C3_prev"] = findAtom(config, "C3", resNum - 1, i, "backward")
				template2 = copy.deepcopy(template)
				for key in template:
					template2[key].pos, _ = config.interatom_vector(template[keyAtom].num, template[key].num)
				template = template2
				A_pos = vector3d.Vector3d()
				A_mass = 0.
				for key in A_lists[unitType]:
					dA = copy.deepcopy(template[key].pos)
					dA.scale(A_lists[unitType][key])
					A_pos += dA
					A_mass += A_lists[unitType][key]
				A_pos.scale(1./A_mass)
				A = Atom()
				A.pos = A_pos
				A.type = "A"
				template["A"] = A
				B_pos = vector3d.Vector3d()
				B_mass = 0.
				for key in B_lists[unitType]:
					dB = copy.deepcopy(template[key].pos)
					dB.scale(B_lists[unitType][key])
					B_pos += dB
					B_mass += B_lists[unitType][key]
				B_pos.scale(1./B_mass)
				B = Atom()
				B.pos = B_pos
				B.type = "B"
				template["B"] = B
				if unitType != "SLS":
					C3_prev_pos = template["C3_prev"].pos
					unitVector = A_pos - C3_prev_pos
				else:
					C3_pos = template["C3"].pos
					unitVector = C3_pos - A_pos
					_, unitVector2, __ = getUnitStructure(config, "SLM", atomNumOffset)
					unitVector += unitVector2
					A_next_pos = A_pos + unitVector
					A_next = Atom()
					A_next.pos = A_next_pos
					template["A_next"] = A_next
				template_counter += 1
				for key in template:
					if key != "O1_next":
						atomNum = template[key].num
						if atomNum != None and atomNum > maxAtomNum:
							maxAtomNum = atomNum
				return template, unitVector, maxAtomNum
#Find key atom, unit type == needed unit type; search for atoms with types from list, res_num = res_num(SLS) + 1; for each atom find image nearest to the C3_prev atom; trim pbc; determine C3 - A distance; add coordinates to database;

def templateToConfig(templateName):
#Read template configuration from file and return it as a Config
	s = open(templateName, 'r').read()
	resType = s.split()[2][:3]
	coords = eval(s)
	config = Config()
	for key in coords:
		atom = Atom()
		atom.is_hetatm = True
		atom.type = re.split('(\d+)', key)[0]
		atom.id = key
		atom.pos = NpArray2Vector3d(coords[key])
		atom.res_type = resType
		config.insert_atom(atom)
	for line in s.splitlines():
		linesplit = line.split()
		if len(linesplit) >= 2 and linesplit[1] == 'bond':
			bond = Bond()
			bond.atom1 = config.atom_by_id(linesplit[2])
			bond.atom2 = config.atom_by_id(linesplit[3])
			config.insert_bond(bond)
	return config
