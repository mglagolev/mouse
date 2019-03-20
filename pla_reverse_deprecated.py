def AlignB(config, coords, B_pos, rotation_center, rotation_vector, n_rot):
	B_discrepancies = []
	rotation_step = 2. * np.pi / n_rot
	rotation_group = { rotation_center : coords[rotation_center], "B" : coords["B"] }
	for i in range(n_rot):
		rotate_group_around(rotation_center, rotation_group, rotation_vector, rotation_step)
		B_discrepancy = vector_pbc_trim(B_pos - NpArray2Vector3d(rotation_group["B"]), config.box()).length()
		B_discrepancies.append([B_discrepancy, i])
	B_discrepancies.sort()
	angle = rotation_step * B_discrepancies[0][1]
	rotate_group_around(rotation_center, coords, rotation_vector, angle)
	return angle

def scale_along_vector(key0, scaling_vector, coords, factor):
	coords_key0 = coords[key0]
	sc_vec_norm = scaling_vector / np.linalg.norm(scaling_vector)
	scaling_vector_sph = Cart2Sph(scaling_vector)
	rotation_vector = np.cross([0., 0., 1.], scaling_vector)
	for key in coords:
		coords[key] = coords_key0 + (np.dot(rotation_matrix(rotation_vector, np.pi - scaling_vector_sph[1]), coords[key] - coords_key0))
		coords[key][2] = coords_key0[2] + (coords[key][2] - coords_key0[2]) * factor
		coords[key] = coords_key0 + (np.dot(rotation_matrix(rotation_vector, scaling_vector_sph[1] - np.pi), coords[key] - coords_key0))


#Rotate around O1-A so that B is in the XY plane
def rotate_B_into_XY(coords):
	rotation_vector = coords["A"] - coords["O1"]
	rotation_angle = math.pi - np.arctan((coords["B"][2] - coords["A"][2])/ (coords["B"][1] - coords["A"][1]))
	for key in coords:
		coords[key] = coords["O1"] + (np.dot(rotation_matrix(rotation_vector,rotation_angle), coords[key] - coords["O1"]))
	return coords


def RestoreByAA(config, coords, A1, B1, A2):
	_coords = coords.copy()
	#This function should return a dictionary containing atom names in AA representation as key and Atom objects as values
	#2. Determine polar coodinates of A1 A2 vector
	sph = Cart2Sph(Vector3d2NpArray(vector_pbc_trim(A2.pos - A1.pos, config.box())))
	#print "CG A-A_next: ", vector_pbc_trim(A2.pos - A1.pos, config.box()), sph
	#3. Rotate coords so that A1 A2 is parallel to coarse-grained representation
	_coords = rotate_group_around("A", _coords, np.array([0.,1.,0.]), -1. * sph[1])
	_coords = rotate_group_around("A", _coords, np.array([0.,0.,1.]), 1. * sph[2])
	#print	"After phi rotation:", "A: ", _coords["A"], " A_next: ", _coords["A_next"], Cart2Sph(_coords["A_next"] - _coords["A"])
	#4. Scale coords so that A2(coords) coincides with coarse-grained A2
	CG_A_A_next = sph[0]
	#CG_A_A_next = np.linalg.norm(Vector3d2NpArray(A2.pos - A1.pos))
	AA_A_A_next = np.linalg.norm(_coords["A_next"] - _coords["A"])
	#AA_A_A_next = sph[0]
	_coords = scale_group(_coords, CG_A_A_next / AA_A_A_next)
	print "Scaling: ", CG_A_A_next / AA_A_A_next
	#1. Translate the coordinates by A1 position vector
	_coords = translate_group(_coords, Vector3d2NpArray(A1.pos))
	#5. Rotate coords around A1-A2 with a small step. Determine rotation which corresponds to minimum distance between B(coords) and coarse-grained B. Rotate coords by this angle
	rotation_vector = _coords["A_next"] - _coords["A"]
	n_rot = 360
	_coords = AlignB(config, _coords, B1, "A", rotation_vector, n_rot)
	return coords

def RestoreByA(config, coords, O1, A, B):
	_coords = coords.copy()
	#This function should return a dictionary containing atom names in AA representation as key and Atom objects as values
	#2. Determine polar coordinates of O1 A vector
	sph = Cart2Sph(Vector3d2NpArray(vector_pbc_trim(A.pos - O1.pos, config.box())))
	#print "CG O1-A: ", vector_pbc_trim(A.pos - O1.pos, config.box()), sph
	#3. Rotate coords so that O1-A is parallel to coarse-grained representation
	_coords = rotate_group_around("O1", _coords, np.array([0.,1.,0.]), -1. * sph[1])
	_coords = rotate_group_around("O1", _coords, np.array([0.,0.,1.]), 1. * sph[2])
	#print	"After phi rotation:", "O1: ", _coords["O1"], " A: ", _coords["A"], Cart2Sph(_coords["A"] - _coords["O1"])
	#4. Scale coords so that A(coords) is parallel to coarse-grained A
	CG_O1_A = sph[0]
	#print "Inside: SPH[0]: ", sph[0]
	#print "O1: ", O1.pos
	#print "A: ", A.pos
	#CG_O1_A = np.linalg.norm(Vector3d2NpArray(A.pos - O1.pos))
	AA_O1_A = np.linalg.norm(_coords["A"] - _coords["O1"])
	#print "Inside: SLM Coords O1-A: ", _coords["A"] - _coords["O1"], "Length: ", np.linalg.norm(_coords["A"] - _coords["O1"])
	#AA_O1_A = sph[0]
	_coords = scale_group(_coords, CG_O1_A / AA_O1_A)
	print "Scaling: ", CG_O1_A / AA_O1_A
	#1. Add coordinates of O1 atom to all atoms in the coords.
	_coords = translate_group(_coords, Vector3d2NpArray(O1.pos))
	#5. Rotate coords around O1-A with a small step. Determine rotation which corresponds to minimum distance between B(coords) and coarse-grained B. Rotate coords by this angle
	rotation_vector = _coords["A"] - _coords["O1"]
	n_rot = 360
	_coords = AlignB(config, _coords, B, "O1", rotation_vector, n_rot)
	return _coords
