import lammpack_types
import vector3d

def AtomFromPdbLine(line):
	"""Returns an Atom object from an atom line in a pdb file."""
	atom = lammpack_types.Atom()
	if line.startswith('HETATM'):
		atom.is_hetatm = True
	else:
		atom.is_hetatm = False
	atom.num = int(line[6:11])
	atom.type = line[12:16].strip(" ")
	element = ''
	for c in line[12:16]:
		if not c.isdigit() and c != " ":
			element += c
	atom.element = element[0]
	atom.res_type = line[17:20]
	atom.mol_id = line[21]
	atom.res_num = int(line[22:26])
	atom.res_insert = line[26]
	if atom.res_insert == " ":
		atom.res_insert = ""
	x = float(line[30:38])
	y = float(line[38:46])
	z = float(line[46:54])
	atom.pos.set(x, y, z)
	try:
		atom.occupancy = float(line[54:60])
	except:
		atom.occupancy = 100.0
	try:
		atom.bfactor = float(line[60:66])
	except:
		atom.bfactor = 0.0
	return atom

def BondsFromPdbLine(line, config):
	"""Returns bond list from CONECT line in a pdb file."""
	bonds = []
	bond_num = config.n_bond()
	bond_atom1 = int(line[6:11])
	bond_atom2_list = []
	str_pos = 11
	while True:
		try:
			bond_atom2 = int(line[str_pos:str_pos+5])
			str_pos += 5
			bond_atom2_list.append(bond_atom2)
		except ValueError:
			break
	for bond_atom2 in bond_atom2_list:
		bond = lammpack_types.Bond()
		bond_num += 1
		bond.num = bond_num
		bond.atom1.num = bond_atom1
		bond.atom2.num = bond_atom2
		try:
			type1 = config.atom_by_num(bond_atom1).type
			type2 = config.atom_by_num(bond_atom2).type
			types = [type1, type2]
			types.sort()
			bond.type = str(types[0]) + str(types[1])
		except: pass
		bonds.append(bond)
	return bonds

def BoxFromPdbLine(line):
	_box = vector3d.Vector3d()
	linesplit = line.split()
	_box.x = float(linesplit[1])
	_box.y = float(linesplit[2])
	_box.z = float(linesplit[3])
	return _box
