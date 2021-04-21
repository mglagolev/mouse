import lammpack_types
import vector3d

def AtomFromGroLine(line):
	"""Returns an Atom object from an atom line in a .gro file."""
	atom = lammpack_types.Atom()
	atom.res_num = int(line[0:5])
	atom.mol_id = atom.res_num	#Residue number used as chain number
	atom.res_type = line[5:10].strip()
	atom.type = line[10:15].strip()
	atom.num = int(line[15:20])
	atom.pos.x = float(line[20:28])
	atom.pos.y = float(line[28:36])
	atom.pos.z = float(line[36:44])
	try:
		atom.vel.x = float(line[44:52])
		atom.vel.y = float(line[52:60])
		atom.vel.z = float(line[60:68])
	except: pass
	return atom
