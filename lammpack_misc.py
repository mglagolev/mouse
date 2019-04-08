#!/usr/bin/env python3

import sys
import lammpack_types
from lammpack_types import Config
import vector3d

def determine_data_type( filename ):
	ftype_extension = 'unknown'
	ftype_content = 'unknown'
#Try to determine file type by extension	
	if len(filename) >= 4:
		ext = filename[-4:].lower()
		if ext == '.ent' or ext == '.pdb':
			ftype_extension = 'pdb'
		elif ext == '.gro':
			ftype_extension = 'gro'
	elif len(filename) >= 7:
		ext = filename[-7:].lower()
		if ext == '.lammps':
			ftype_extension = 'lammps'
#Try to determine file type by keywords
	infile = open( filename, 'rb' )
	header = infile.readline()
	
	if len(header) >= 16:	
		if header[:16].lower() == 'lammps data file' or header[:18].lower() == 'lammps description':
			ftype_content = 'lammps_data'
	if len(header) >= 6:
		h6low = header[:6].lower()
		if h6low == 'model ' or h6low == 'atom  ' or h6low == 'hetatm' or h6low == 'conect':
			ftype_content = 'pdb'
	if len(header) >= 5:
		if header[:5].lower() == 'item:':
			ftype_content = 'lammps_dump'
	if ftype_extension != 'unknown':
		return ftype_extension
	else:
		return ftype_content


def read_data_typeselect( filename, ftype=None, options = {} ):
	"""
	Read data file according to its file type

	Type-specific options can be passed to the file-reading function through "options" key-value pairs
	"""
	if ftype is None or ftype == 'auto':
		ftype = determine_data_type( filename )
	readfunc = {'lammps_dump':1, 'lammps_data':2, 'pdb':3, 'gro':4}
	if ftype in readfunc:
		try: options = options[ftype]
		except KeyError: options = {}
		config = Config()
		if ftype == 'lammps_dump':
			config.read_lmp_dump(filename, **options)
		elif ftype == 'lammps_data':
			config.read_lmp_data(filename, **options)
		elif ftype == 'pdb':
			config.read_pdb(filename, **options)
		elif ftype == 'gro':
			config.read_gro(filename, **options)
		return config
	else:
		error_str = "No handler function for file type " + ftype
		raise NameError(error_str)

def vector_pbc_trim( vector, box ):
#Calculate the values of the vector between the nearest replicas taking into account the periodic boundary conditions
	x = vector.x - box.x * int(round(vector.x / box.x))
	y = vector.y - box.y * int(round(vector.y / box.y))
	z = vector.z - box.z * int(round(vector.z / box.z))
	_vector_pbc_trim = vector3d.Vector3d()
	_vector_pbc_trim.set( x, y, z )
	return _vector_pbc_trim

def selectConstraintsByAtom(config, atomsList):
	newConstraints = { "bonds" : [], "angles" : [], "dihedrals" : [] }
	atomNumbersList = [ x.num for x in atomsList ]
	for bond in config.bonds():
		try:
			atomNumbersList.index(bond.atom1.num)
			atomNumbersList.index(bond.atom2.num)
			newConstraints["bonds"].append(bond)
		except ValueError: pass
	for angle in config.angles():
		try:
			atomNumbersList.index(angle.atom1.num)
			atomNumbersList.index(angle.atom2.num)
			atomNumbersList.index(angle.atom3.num)
			newConstraints["angles"].append(angle)
		except ValueError: pass
	for dihedral in config.dihedrals():
		try:
			atomNumbersList.index(dihedral.atom1.num)
			atomNumbersList.index(dihedral.atom2.num)
			atomNumbersList.index(dihedral.atom3.num)
			atomNumbersList.index(dihedral.atom4.num)
			newConstraints["dihedrals"].append(dihedral)
		except ValueError: pass
	return newConstraints


def selectByAtomFields(config, selectionCriteria):
	"""
	Create new configuration with atoms selected by fields matching the criteria.

	selectionCriteria: { fieldName1 : [ value1, value2, ... ], fieldName2 : [ value1, value2, ...] }

	The bonds, angles and dihedrals with all atoms satisfying the criteria are also included to the new configuration.
	"""
	newConfig = Config()
	newConfig.set_box(config.box())
	newConfig.set_box_center(config.box_center())
	for atom in config.atoms():
		try:
			for fieldName in selectionCriteria:
				fieldValue = getattr(atom, fieldName)
				selectionCriteria[fieldName].index(fieldValue)
			newConfig.insert_atom(atom)
		except ValueError: pass
	newConstraints = selectConstraintsByAtom(config, newConfig.atoms())
	newBonds, newAngles, newDihedrals = newConstraints["bonds"], newConstraints["angles"], newConstraints["dihedrals"]
	for bond in newBonds: newConfig.insert_bond(bond)
	for angle in newAngles: newConfig.insert_angle(angle)
	for dihedral in newDihedrals: newConfig.insert_dihedral(dihedral)
	return newConfig	 


def select_chains(config, chainList):
	return selectByAtomFields(config, { "mol_id" : chainList })


def select_rectangular(config, cutx = False, xrelmin = 0., xrelmax = 1., cuty = False, yrelmin = 0., yrelmax = 1., cutz = False, zrelmin = 0., zrelmax = 1.):
	"""
	Create new configuration with atoms contained in the rectangular subcell. The coordinates are specified relative to the original cell.

	The bonds, angles and dihedrals with all atoms contained in the new subcell are also included to the new configuration.
	"""
	newConfig = Config()
	newConfig.set_box(config.box())
	newConfig.set_box_center(config.box_center())
	for atom in config.atoms():
		xrel = atom.pos.x / config.box().x
		yrel = atom.pos.y / config.box().y
		zrel = atom.pos.z / config.box().z
		if ( cutx == False or (xrel >= xrelmin and xrel <= xrelmax)) and ( cuty == False or (yrel >= yrelmin and yrel <= yrelmax)) and (cutz == False or (zrel >= zrelmin and zrel <= zrelmax)):
			newConfig.insert_atom(atom)
	newConstraints = selectConstraintsByAtom(config, newConfig.atoms())
	newBonds, newAngles, newDihedrals = newConstraints["bonds"], newConstraints["angles"], newConstraints["dihedrals"]
	for bond in newBonds: newConfig.insert_bond(bond)
	for angle in newAngles: newConfig.insert_angle(angle)
	for dihedral in newDihedrals: newConfig.insert_dihedral(dihedral)
	return newConfig


def GetMoleculeAtomlists(config, atom_types=["All"]):
	natom = config.n_atom()
	molecules_index = []
	molecules_atoms = []
	for i in range(natom):
		atom = config.atom(i)
		mol = atom.mol_id
		check_types = True
		try:
			if atom_types[0] == "All":
				check_types = False
			else:
				atom_type = atom.type
		except IndexError: atom_type = atom.type
		if check_types:
			try:
				atom_types.index(atom_type)
			except ValueError:
				continue
		try:
			molecules_atoms[molecules_index.index(mol)].append(atom.num)
		except ValueError:
			molecules_index.append(mol)
			molecules_atoms.append([atom.num])
	return molecules_index, molecules_atoms


