import lammpack_types

def BondFromLmpDataLine(line):
	bond = lammpack_types.Bond()
	linesplit = line.split()
	bond.num = int(linesplit[0])
	bond.type = linesplit[1]
	atom1 = int(linesplit[2])
	atom2 = int(linesplit[3])
	if atom1 < atom2:
		bond.atom1.num = atom1
		bond.atom2.num = atom2
	else:
		bond.atom1.num = atom2
		bond.atom2.num = atom1
	return bond

def AtomFromLmpDataLine(line):
	atom = lammpack_types.Atom()
	linesplit = line.split()
	atom.num = int(linesplit[0])
	atom.mol_id = linesplit[1]
	atom.type = linesplit[2]
	#atom.res_num = atom.num
	atom.res_num = int(atom.mol_id)
	atom.pos.x = float(linesplit[3])
	atom.pos.y = float(linesplit[4])
	atom.pos.z = float(linesplit[5])
	atom.pbc.x = int(linesplit[6])
	atom.pbc.y = int(linesplit[7])
	atom.pbc.z = int(linesplit[8])
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
