#!/usr/bin/env python3

"""
Routines to read/write to and from various LAMMPS dump files
"""

import os
import vector3d
from lammpack_lmp_functions import *
from lammpack_pdb_functions import *
from gro_functions import *
from clustering_functions import *
import functools

class Atom:
	def __init__(self):
		self.is_hetatm = False
		self.pos = vector3d.Vector3d()
		self.vel = vector3d.Vector3d()
		self.mass = 0.0
		self.type = ""
		self.element = ""
		self.id = ""
		self.mol_id = ""
		self.res_type = ""
		self.res_num = 0
		self.res_insert = ""
		self.bfactor = 0.0
		self.occupancy = 0.0
		self.num = 0
		self.pbc = vector3d.Vector3d(0,0,0)
		self.neighbors = []
		self.inum = 0 #Internal number within a frame

	def pdb_str(self):
		if self.is_hetatm:
			field = "HETATM"
		else:
			field = "ATOM  "
		return "%6s%5s %4s %4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f         %2s%2s" \
            % (field, self.num, 
               pad_atom_type(self.type),
               self.res_type, self.mol_id[:1],
               self.res_num, self.res_insert,
               self.pos.x, self.pos.y, self.pos.z,
               self.occupancy, self.bfactor, "  ", "  ")
               
	def __str__(self):
		return "%s%s-%s (% .1f % .1f % .1f)" \
            %  (self.res_type, self.res_num, 
                self.type, self.pos.x, 
                self.pos.y, self.pos.z)

	def gro_str(self):
		if self.res_type == '':
			self.res_type = 'DUMMY'
		return "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" \
		% ( self.res_num, self.res_type[:5],
		self.type[:5], self.num,
		self.pos.x, self.pos.y, self.pos.z,
		self.vel.x, self.vel.y, self.vel.z)

	def add_neighbor(self, atom):
		self.neighbors.append(atom)

	def remove_neighbor(self, atom):
		self.neighbors.remove(atom)


class Bond:
	def __init__(self):
		self.atom1 = Atom()
		self.atom2 = Atom()
		self.type = ""
		self.dist = 0.0
		self.energy = 0.0
		self.num = 0

	def pdb_str(self):
		field = "CONECT"
		return "%6s%5s%5s" % (field, self.atom1.num, self.atom2.num)


class Angle:
	def __init__(self):
		self.atom1 = Atom()
		self.atom2 = Atom()
		self.atom3 = Atom()
		self.num = 0
		self.type = ""


class Dihedral:
	def __init__(self):
		self.atom1 = Atom()
		self.atom2 = Atom()
		self.atom3 = Atom()
		self.atom4 = Atom()
		self.num = 0
		self.type = ""

class Molecule:
	def __init__(self):
		self.atoms = []
		self.bonds = []
		self.angles = []
		self.dihedrals = []

class Config:
	def __init__(self):
		self._atoms = []
		self._bonds = []
		self._angles = []
		self._dihedrals = []
		self._molecules = []
		self._box = vector3d.Vector3d()
		self._box_center = vector3d.Vector3d()
		self.title = ""

	def box(self):
		return self._box

	def set_box(self, box):
		self._box = box

	def box_center(self):
		return self._box_center

	def set_box_center(self, box_center):
		self._box_center = box_center

	def n_atom(self):
		return len(self._atoms)

	def n_bond(self):
		return len(self._bonds)

	def atoms(self):
		return self._atoms

	def atom(self, i):
		return self._atoms[i]

	def atom_by_num(self, i):
		for atom in self._atoms:
			if atom.num == i:
				return atom
	def atom_by_inum(self, i):
		atomIndex = [ x.inum for x in self._atoms ].index(i)
		return self._atoms[atomIndex]

	def atom_by_id(self, s):
		for atom in  self._atoms:
			if atom.id == s:
				return atom

	def bonds(self):
		return self._bonds

	def bond(self, i):
		return self._bonds[i]

	def bond_by_num(self, i):
		return self._bonds[[ x.num for x in self._bonds ].index(i)]

	def angles(self):
		return self._angles

	def dihedrals(self):
		return self._dihedrals
    
	def clear(self):
		for atom in self._atoms:
			del atom
		del self._atoms[:]

	def transform(self, matrix):
		for atom in self._atoms:
			atom.pos.transform(matrix)

	def convertCoordsCutToUncut(self):
		for atom in self._atoms:
			atom.pos.x = atom.pos.x + self._box.x * atom.pbc.x
			atom.pos.y = atom.pos.y + self._box.y * atom.pbc.y
			atom.pos.z = atom.pos.z + self._box.z * atom.pbc.z
			

	def insert_atom(self, atom):
		self._atoms.append(atom)
    
	def erase_atom(self, atom_type):
		for atom in self._atoms:
			if atom.type == atom_type:
				self._atoms.remove(atom)
				del atom
				return
	def insert_bond(self, bond):
		self._bonds.append(bond)

	def insert_angle(self, angle):
		self._angles.append(angle)

	def insert_dihedral(self, dihedral):
		self._dihedrals.append(dihedral)

	def read_pdb(self, fname, assignMolecules = { "type" : "chainId", "cluster" : False }):
		"""
		Create new configuration from .pdb file.

		assignMolecules:
			"type":
				"chainId": the molecules will be determined from "chainID" field (position 22 in "ATOM/HETATM" record).
				"resSeq": the molecules will be determined from "resSeq" field (positions 23-26 in "ATOM/HETATM" record).
			"cluster":
				False: If the corresponding field is empty, "mol_id" will be left empty.
				True: In case the field is empty for at least one atom, all "mol_id"s will be assigned by cluster determination algorithm based on atom connectivity.
		"""
		self.clear()
		inum = 1
		for line in open(fname, 'r').readlines():
			if line.startswith("ATOM") or line.startswith("HETATM"):
				atom = AtomFromPdbLine(line, inum, assignMolecules)
				if len(self._atoms) == 1:
					self.id = atom.mol_id
				self.insert_atom(atom)
				inum += 1
			if line.startswith("CONECT"):
				bonds = BondsFromPdbLine(line, self)
				for bond in bonds:
					self.insert_bond(bond)
			if line.startswith("ENDMDL"):
				continue
			if line.startswith("CRYST1"):
				box = BoxFromPdbLine(line)
				self.set_box(box)
			if line.startswith("TITLE"):
				self.title = line[10:].rstrip()
		if assignMolecules["cluster"]:
			for atom in self._atoms:
				if atom.mol_id == "" or atom.mol_id == " ":
					assignMoleculesFromBonds(self)
					break

	def write_pdb(self, pdb, hide_pbc_bonds = False):
		f = open(pdb, 'w')
		n_atom = 0
		f.write("TITLE     " + self.title + "\n")
		f.write("%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f\n" \
            % ("CRYST1", self.box().x, self.box().y, self.box().z, 90., 90., 90.))
		f.write("MODEL     %4i\n" % (1))
		for atom in sorted(self._atoms, cmp=cmp_atom):
			n_atom += 1
			atom.num = n_atom
			f.write(atom.pdb_str() + '\n')
		f.write("TER\n")
		f.write("ENDMDL\n")
		for bond in sorted(self._bonds, cmp=cmp_bond):
			visible= True
			if hide_pbc_bonds:
				atom1 = self.atom_by_num(bond.atom1.num)
				atom2 = self.atom_by_num(bond.atom2.num)
				box_x, box_y, box_z = self.box().x, self.box().y, self.box().z
				dx, dy, dz = abs(atom2.pos.x - atom1.pos.x), abs(atom2.pos.y - atom1.pos.y), abs(atom2.pos.z - atom1.pos.z)
				if (box_x > 0 and dx / box_x > 0.5) or (box_y > 0 and dy / box_y > 0.5) or (box_z > 0 and dz / box_z > 0.5):
					visible = False
			if visible:
				f.write(bond.pdb_str() + '\n')
		f.close()

	def read_lmp_dump(self, lmp_dump):
		raise NameError('Reading LAMMPS dump is not yet implemented')

	def read_lmp_data(self, lmp_data):
		f = open(lmp_data, 'r')
		while True:
			line=f.readline()
			if len(line) == 0: break
			linesplit = line.split()
			if len(linesplit) == 0:
				continue
			else:
#Read bonds
				if linesplit[0] == "Bonds":
					line = f.readline()
					for i in xrange(nbond):
						bond = BondFromLmpDataLine( f.readline() )
						self.insert_bond(bond)

#Read angles
				elif linesplit[0] == "Angles":
					line = f.readline()
					for i in xrange(nangle):
						angle = AngleFromLmpDataLine( f.readline() )
						self.insert_angle(angle)

#Read atoms	
				elif linesplit[0] == "Atoms":
					line = f.readline()
					for i in xrange(natom):
						atom = AtomFromLmpDataLine( f.readline() )
						self.insert_atom(atom)

#Read dihedrals
				elif linesplit[0] == "Dihedrals":
					line = f.readline()
					for i in xrange(ndihedrals):
						dihedral = DihedralFromLmpDataLine( f.readline() )
						self.insert_dihedral(dihedral)

				elif len(linesplit) >= 4 and linesplit[2] == 'xlo' and linesplit[3] == 'xhi':
					xlo, xhi = float(linesplit[0]), float(linesplit[1])
					self._box.x = abs(xhi - xlo)
					self._box_center.x = (xhi + xlo) / 2.
				elif len(linesplit) >= 4 and linesplit[2] == 'ylo' and linesplit[3] == 'yhi':
					ylo, yhi = float(linesplit[0]), float(linesplit[1])
					self._box.y = abs(yhi - ylo)
					self._box_center.y = (yhi + ylo) / 2.
				elif len(linesplit) >= 4 and linesplit[2] == 'zlo' and linesplit[3] == 'zhi':
					zlo, zhi = float(linesplit[0]), float(linesplit[1])
					self._box.z = abs(zhi - zlo)
					self._box_center.z = (zhi + zlo) / 2.
#Read header: n_bonds and n_angles
				elif len(linesplit) >=2 and linesplit[1] == "atoms":
					natom = int(linesplit[0])
				elif len(linesplit) >=2 and linesplit[1] == "bonds":
					nbond = int(linesplit[0])
				elif len(linesplit) >=2 and linesplit[1] == "angles":
					nangle = int(linesplit[0])
				elif len(linesplit) >=2 and linesplit[1] == "dihedrals":
					ndihedrals = int(linesplit[0])
				else:
					if linesplit[0] == "LAMMPS":
						continue
		f.close()

	def bond_vector_by_num(self, bond_num):
		bond = next((x for x in self._bonds if x.num == bond_num), None)
		return self.interatom_vector(bond.atom1.num, bond.atom2.num)

	def interatom_vector(self, atom1_num, atom2_num):
		atom1 = next((x for x in self._atoms if x.num == atom1_num), None)
		atom2 = next((x for x in self._atoms if x.num == atom2_num), None)
		x1, y1, z1 = atom1.pos.x, atom1.pos.y, atom1.pos.z
		x2, y2, z2 = atom2.pos.x, atom2.pos.y, atom2.pos.z
		x = x2 - x1
		y = y2 - y1
		z = z2 - z1
		xcenter = (x2 + x1) / 2.
		ycenter = (y2 + y1) / 2.
		zcenter = (z2 + z1) / 2.
		if self._box.x > 0:
			x_trimmed = x - self._box.x * int(round(x / self._box.x))
			xcenter_trimmed = xcenter -  0.5 * self._box.x * int(round(x / self._box.x))
			xcenter_trimmed = xcenter_trimmed - self._box.x * int(round(xcenter_trimmed / self._box.x))
		else:
			x_trimmed = x
			xcenter_trimmed = xcenter
			xcenter_trimmed = xcenter_trimmed - self._box.x * int(round(xcenter_trimmed / self._box.x))
		if self._box.y > 0:
			y_trimmed = y - self._box.y * int(round(y / self._box.y))
			ycenter_trimmed = ycenter -  0.5 * self._box.y * int(round(y / self._box.y))
			ycenter_trimmed = ycenter_trimmed - self._box.y * int(round(ycenter_trimmed / self._box.y))
		else:
			y_trimmed = y
			ycenter_trimmed = ycenter
			ycenter_trimmed = ycenter_trimmed - self._box.y * int(round(ycenter_trimmed / self._box.y))
		if self._box.z > 0:
			z_trimmed = z - self._box.z * int(round(z / self._box.z))
			zcenter_trimmed = zcenter -  0.5 * self._box.z * int(round(z / self._box.z))
			zcenter_trimmed = zcenter_trimmed - self._box.z * int(round(zcenter_trimmed / self._box.z))
		else:
			z_trimmed = z
			zcenter_trimmed = zcenter
			zcenter_trimmed = zcenter_trimmed - self._box.z * int(round(zcenter_trimmed / self._box.z))
		_bond_vector = vector3d.Vector3d()
		_bond_vector.set(x_trimmed, y_trimmed, z_trimmed)
		_bondcenter_vector = vector3d.Vector3d()
		_bondcenter_vector.set(xcenter_trimmed, ycenter_trimmed, zcenter_trimmed)
		return _bond_vector, _bondcenter_vector

	def read_gro(self, fname):
		self.clear()
		f = open(fname, 'r')
		self.title = f.readline().strip()
		natoms = int(f.readline())
		for i in range(natoms):
			atom = AtomFromGroLine(f.readline())
			self.insert_atom(atom)
		box_array = map(float, f.readline().split())
		_box = vector3d.Vector3d(box_array[0], box_array[1], box_array[2])
		self.set_box(_box)
		f.close()
			
	def write_gro(self, fname):			
		f = open(fname, 'w')
		f.write(self.title + "\n")
		f.write(str(self.n_atom()) + "\n")
		for atom in sorted(self._atoms, cmp=cmp_atom):
			f.write(atom.gro_str() + '\n')
		f.write(" " + str(self._box.x).rjust(9) + " " + str(self._box.y).rjust(9) + " " + str(self._box.z).rjust(9))
		f.close()

	def trim(self):
		sorted(self._atoms, key = functools.cmp_to_key(cmp_atom))
		for bond in self._bonds:
			if bond.atom1.num > bond.atom2.num:
				atom1, atom2 = bond.atom1.num, bond.atom2.num
				bond.atom1.num, bond.atom2.num = atom2, atom1
		sorted(self._bonds, key = functools.cmp_to_key(cmp_bond))

	def reassign_inums(self):
		inum = 1
		for atom in self._atoms:
			atom.inum = inum
			inum += 1 

def cmp_atom(a1, a2):
	if a1.num < a2.num:
		return -1
	else:
		return 0

def cmp_bond(b1, b2):
	if b1.atom1.num < b2.atom1.num:
		return -1
	else:
		return 0

module_dir = os.path.dirname(__file__)
radii_fname = os.path.join(module_dir, "radii.txt")
f = open(radii_fname, 'r')
radii = eval(f.read())
f.close()
two_char_elements = [el for el, r in radii.items() if len(el) == 2]

def pad_atom_type(in_atom_type):
  atom_type = in_atom_type
  if len(atom_type) == 1:
    atom_type = " %s  " % atom_type
  elif len(atom_type) == 2:
    atom_type = " %s " % atom_type
  elif len(atom_type) == 3:
    if atom_type[0].isdigit():
      atom_type = "%s " % atom_type
    else:
      atom_type = " %s" % atom_type
  return atom_type
