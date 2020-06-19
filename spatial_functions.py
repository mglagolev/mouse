#!/usr/bin/env python

from lammpack_misc import *
from lammpack_types import *
from vector3d import *
import numpy as np
import sys

def CalculateEnd2EndDistance(config, atom_types):
	molecules_index, molecules_atoms = GetMoleculeAtomlists(config, atom_types)
	average_end2end_dist = 0.
	end2end_dist_counter = 0
	for i in range(len(molecules_index)):
		molecule_id = molecules_index[i]
		molecule_atoms = molecules_atoms[i]
		firstatom = config.atom_by_num(min(molecule_atoms))
		lastatom = config.atom_by_num(max(molecule_atoms))
		end2end_vector_untrimmed = Vector3d(lastatom.pos.x - firstatom.pos.x, lastatom.pos.y - firstatom.pos.y, lastatom.pos.z - firstatom.pos.z)
		end2end_vector = vector_pbc_trim(e2e_vector_untrimmed, config.box())
		end2end_dist = e2e_vector.length()
		average_end2end_dist += end2end_dist
		end2end_dist_counter += 1
	average_end2end_dist /= end2end_dist_counter
	return average_end2end_dist		

def CalculateLoopPitch(config, atom1_number, atom2a_number, atom2b_number, period):
	period_trunc = int(period)
	period_fraction = period - period_trunc
	atom1 = config.atom_by_num(atom1_number)
	atom2a = config.atom_by_num(atom2a_number)
	atom2b = config.atom_by_num(atom2b_number)
	x1, y1, z1 = atom1.pos.x, atom1.pos.y, atom1.pos.z
	x2a, y2a, z2a = atom2a.pos.x, atom2a.pos.y, atom2a.pos.z
	x2b, y2b, z2b = atom2b.pos.x, atom2b.pos.y, atom2b.pos.z
	x2 = x2a * (1. - period_fraction) + x2b * period_fraction
	y2 = y2a * (1. - period_fraction) + y2b * period_fraction
	z2 = z2a * (1. - period_fraction) + z2b * period_fraction
	pitch_vector_untrimmed = Vector3d(x2 - x1, y2 - y1, z2 - z1)
	pitch_vector = vector_pbc_trim(pitch_vector_untrimmed, config.box())
	pitch = pitch_vector.length()
	return pitch


def CalculateHelixPitch(config, atom_types, period):
	period_trunc = int(period)
	period_fraction = period - period_trunc
	molecules_index, molecules_atoms = GetMoleculeAtomlists(config, atom_types)
	molecules_average_pitch = 0.
	molecule_counter = 0
	for i in range(len(molecules_index)):
		molecule_id = molecules_index[i]
		molecule_atoms = molecules_atoms[i]
		molecule_atoms.sort()
		average_pitch = 0.
		pitch_counter = 0
		for j in range(len(molecule_atoms) - period_trunc - 1):
			pitch = CalculateLoopPitch(config, molecule_atoms[j], molecule_atoms[j+period_trunc], molecule_atoms[j+period_trunc+1], period)
			average_pitch += pitch
			pitch_counter += 1
		average_pitch /= pitch_counter
		molecules_average_pitch += average_pitch
		molecule_counter += 1
	molecules_average_pitch /= molecule_counter
	return molecules_average_pitch


def CalculateHelixTubeRadius(config, atom_types, period):
	period_trunc = int(period)
	period_fraction = period - period_trunc
	molecules_index, molecules_atoms = GetMoleculeAtomlists(config, atom_types)
	molecules_r_tube = 0.
	molecules_r_tube_counter = 0
	for i in range(len(molecules_index)):
		molecule_id = molecules_index[i]
		molecule_atoms = molecules_atoms[i]
		molecule_atoms.sort()
		xstart, ystart, zstart = 0., 0., 0.
		for j in range(period_trunc):
			atom = config.atom_by_num(molecule_atoms[j])
			xstart += atom.pos.x
			ystart += atom.pos.y
			zstart += atom.pos.z
		atom = config.atom_by_num(molecule_atoms[period_trunc])
		xstart += atom.pos.x * period_fraction
		ystart += atom.pos.y * period_fraction
		zstart += atom.pos.z * period_fraction
		xstart /= period
		ystart /= period
		zstart /= period
		xend, yend, zend = 0., 0., 0.
		for j in range(period_trunc):
			atom = config.atom_by_num(molecule_atoms[-j - 1])
			xend += atom.pos.x
			yend += atom.pos.y
			zend += atom.pos.z
		atom = config.atom_by_num(molecule_atoms[-1 * period_trunc - 1])
		xend += atom.pos.x * period_fraction
		yend += atom.pos.y * period_fraction
		zend += atom.pos.z * period_fraction
		xend /= period
		yend /= period
		zend /= period
		axis_start = np.array([xstart, ystart, zstart])
		axis_end = np.array([xend, yend, zend])
		r_tube = 0.
		r_tube_counter = 0
		for j in range(period_trunc+1, len(molecule_atoms) - period_trunc - 1):
			atom = config.atom_by_num(molecule_atoms[j])
			bead_position = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
			r_tube += float(np.linalg.norm(np.cross(axis_end - axis_start, axis_start - bead_position))) / float(np.linalg.norm(axis_end - axis_start))
			sys.stderr.write(str(float(np.linalg.norm(np.cross(axis_end - axis_start, axis_start - bead_position))) / float(np.linalg.norm(axis_end - axis_start))) + '\n')
			r_tube_counter += 1
		r_tube /= r_tube_counter
		molecules_r_tube += r_tube
		molecules_r_tube_counter += 1
	molecules_r_tube /= molecules_r_tube_counter
	return molecules_r_tube


def CalculateHelixTubeRadius2(config, atom_types, period):
	period_trunc = int(period)
	halfperiod = period / 2.
	halfperiod_trunc = int(halfperiod)
	halfperiod_fraction = halfperiod - halfperiod_trunc
	quarterperiod = period / 4.
	quarterperiod_trunc = int(quarterperiod)
	molecules_index, molecules_atoms = GetMoleculeAtomlists(config, atom_types)
	molecules_r_tube = 0.
	molecules_r_tube_counter = 0
	for i in range(len(molecules_index)):
		molecule_id = molecules_index[i]
		molecule_atoms = molecules_atoms[i]
		molecule_atoms.sort()
		r_tube = 0.
		r_tube_counter = 0
		for j in range(quarterperiod_trunc, len(molecule_atoms) + quarterperiod_trunc - period_trunc - 1):
			#Calculate pitch value for helix loop starting with atom i
			j_pitch = j - quarterperiod_trunc
			pitch = CalculateLoopPitch(config, molecule_atoms[j_pitch], molecule_atoms[j_pitch+period_trunc], molecule_atoms[j_pitch+period_trunc+1], period)
			#Calculate helix tube diameter for helix half-loop starting with atom i
			alpha = 2. * np.pi / period
			atom1 = config.atom_by_num(molecule_atoms[j])
			atom2a = config.atom_by_num(molecule_atoms[j + halfperiod_trunc])
			atom2b = config.atom_by_num(molecule_atoms[j + halfperiod_trunc + 1])
			vector_ab = np.array([atom2b.pos.x - atom2a.pos.x, atom2b.pos.y - atom2a.pos.y, atom2b.pos.z - atom2a.pos.z])
			r_ab = np.linalg.norm(vector_ab)
			r_tube_est = r_ab / 2. / np.sin(alpha / 2.)
			alpha_a = alpha * halfperiod_fraction
			alpha_b = alpha * (1. - halfperiod_fraction)
			beta_a = (alpha - alpha_a) / 2.
			a = 2. * r_tube_est * np.sin(alpha_a / 2. )
			b = 2. * r_tube_est * np.sin(alpha_b / 2. )
			cos_beta_a = ( a**2 + r_ab**2 - b**2 ) / 2. / a / r_ab
			sin_beta_a = ( 1 - cos_beta_a**2 )**0.5
			upper2 = np.sin(alpha_b/2.)*np.sin(alpha_a/2.)
			lower2 = r_tube_est * np.sin(alpha_a)
			upper = np.cos(alpha_b/2.) - upper2 / lower2
			lower = 1. - a**2 * (np.sin(alpha_b/2.)/r_tube_est/np.sin(alpha_a))**2
			tangential_amendment = a * upper / lower		
			radial_amendment = tangential_amendment * a * np.sin(alpha_b / 2.) / r_tube_est / np.sin(alpha_a)
			vector_1_2a = np.array([atom2a.pos.x - atom1.pos.x, atom2a.pos.y - atom1.pos.y, atom2a.pos.z - atom1.pos.z])
			vector_skewed = vector_1_2a + vector_ab * tangential_amendment / r_ab
			vector_skewed = vector_skewed * ( 1. + radial_amendment / np.linalg.norm(vector_skewed))
			d_tube_sq = (np.linalg.norm(vector_skewed))**2 - (pitch / 2.)**2
			if d_tube_sq >= 0.:
				r_tube += d_tube_sq**0.5 / 2.
				#sys.stderr.write(str(d_tube_sq**0.5 / 2.) + '\n')
				r_tube_counter += 1
		r_tube /= r_tube_counter
		molecules_r_tube += r_tube
		molecules_r_tube_counter += 1
	molecules_r_tube /= molecules_r_tube_counter
	return molecules_r_tube
	
def CalculatePhaseVolume(config, trials, rprobe = -1, error_estimate_runs = 10):
	""" Calculate relative phase volume for each atom type by Monte-Carlo """
	import random
	import ordering_functions as of
	if rprobe == -1: rprobe = config.box().x / float(trials)**(1./3.) / 2.
	types = {}
	n_estimate = int(trials / error_estimate_runs)
	for atom in config.atoms():
		if atom.type not in types:
			types[atom.type] = 0
	npAtomsByType = {}
	for atype in types:
		npAtomsByType[atype] = of.createNumpyAtomsArrayFromConfig(config, allowedAtomTypes = [atype])
	phase = types.copy()
	error_estimate = {}
	stdev = {}
	for atype in types:
		error_estimate[atype] = []
	for itry in range(trials):
		x, y, z = random.random() * config.box().x, random.random() * config.box().y, random.random() * config.box().z
		probe = Atom()
		probe.pos = vector3d.Vector3d(x, y, z)
		neighbors = types.copy()
		for atype in types:
			rdf = of.calculateRdfForReference(config, probe, npAtomsByType[atype], rmin = 0., rmax = rprobe, nbin = 1, excludeSelf = False, sameMolecule = True)
			neighbors[atype] = rdf[0][0]
		total_neighbors = 0
		for atype in neighbors:
			total_neighbors += neighbors[atype]
		if total_neighbors > 0:	
			for atype in neighbors:
				if float(neighbors[atype]) / float(total_neighbors) >= 0.5:
					phase[atype] += 1
					break
		if itry % 100 == 0:
			sys.stderr.write("\rTrials completed: " + str(itry))
		if itry % n_estimate == 0 and itry > 0:
			for atype in types:
				if len(error_estimate[atype]) >= 1:
					error_estimate[atype].append(phase[atype] - np.sum(error_estimate[atype]))
				else:
					error_estimate[atype].append(phase[atype])
	for atype in phase:
		phase[atype] = float(phase[atype]) / float(trials)
		for i in range(len(error_estimate[atype])):
			error_estimate[atype][i] /= n_estimate
		sys.stderr.write(str(error_estimate[atype]) + "\n")
		stdev[atype] = np.std(error_estimate[atype])
	return phase, stdev
