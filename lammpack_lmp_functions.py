import lammpack_types

def BondFromLmpDataLine(line, config):
	bond = lammpack_types.Bond()
	linesplit = line.split()
	bond.num = int(linesplit[0])
	bond.type = linesplit[1]
	atom_nums = [int(linesplit[2]) ,int(linesplit[3])]
	atom_nums.sort()
	bond.atom1 = config.atom_by_num(atom_nums[0])
	bond.atom2 = config.atom_by_num(atom_nums[1])
	bond.num = config.n_bond() + 1
	return bond

def AtomFromLmpDataLine(line):
	atom = lammpack_types.Atom()
	linesplit = line.split()
	nvalues = len(linesplit)
	atom.id = linesplit[0]
	atom.num = int(linesplit[0])
	atom.mol_id = linesplit[1]
	atom.type = linesplit[2]
	#atom.res_num = atom.num
	atom.res_num = int(atom.mol_id)
	if len(linesplit) == 9:
		atom.charge = float(linesplit[3])
	atom.pos.x = float(linesplit[nvalues - 6])
	atom.pos.y = float(linesplit[nvalues - 5])
	atom.pos.z = float(linesplit[nvalues - 4])
	atom.pbc.x = int(linesplit[nvalues - 3])
	atom.pbc.y = int(linesplit[nvalues - 2])
	atom.pbc.z = int(linesplit[nvalues - 1])
	return atom

def AngleFromLmpDataLine(line):
	angle = lammpack_types.Angle()
	linesplit = line.split()
	angle.num = int(linesplit[0])
	angle.type = linesplit[1]
	angle.atom1.num = int(linesplit[2])
	angle.atom2.num = int(linesplit[3])
	angle.atom3.num = int(linesplit[4])
	return angle

def DihedralFromLmpDataLine(line):
	dihedral = lammpack_types.Dihedral()
	linesplit = line.split()
	dihedral.num = int(linesplit[0])
	dihedral.type = linesplit[1]
	dihedral.atom1.num = int(linesplit[2])
	dihedral.atom2.num = int(linesplit[3])
	dihedral.atom3.num = int(linesplit[4])
	dihedral.atom4.num = int(linesplit[5])
	return dihedral
