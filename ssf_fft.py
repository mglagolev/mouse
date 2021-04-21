#Using P. G. Khalatur's SK subroutine
import sys
import math
import cmath
from numpy import zeros
from numpy import fft
import os
import lammpack_types

#-------CORE FUNCTION SK------------------------------------------------
def skfft(config):
#***********************************************************************
#*    static structure factor for bicomponent system
#*
#*    particle positions: x, y, z
#*
#*    The spherical cutoff is used for the reciprocal lattice vectors
#*
#*    rksqmax = ksqmax * twopi**2 / boxmin**2
#*
#***********************************************************************

	box3d = config.box()
	box_center = config.box_center()
	box = min(box3d.x, box3d.y, box3d.z)

	x, y, z = [], [], []

	for atom in config.atoms():
		ux = atom.pos.x
		uy = atom.pos.y
		uz = atom.pos.z

		sx = (ux - box_center.x) / box3d.x + 0.5
		sy = (uy - box_center.y) / box3d.y + 0.5
		sz = (uz - box_center.z) / box3d.z + 0.5

		x.append(sx)
		y.append(sy)
		z.append(sz)

	natom = config.n_atom()

	nmax, kmax, ksqmax = natom, 25, 1024
	pid = str(os.getpid())
	expikr, summ = 0+0j, 0+0j
	histo = zeros((kmax+1, kmax+1, kmax+1), float)
	ns = zeros(ksqmax, float)
	sk = zeros(ksqmax, float)


#.....set initial values
	if 1:
		if not os.path.isfile('sk-' + pid + '.tmp'):
			isk = 0
		else:
			out = open('sk-' + pid + '.tmp', 'r+')
			isk = int(out.readline().split()[0])
 			for i in xrange(ksqmax):
				datum = out.readline().split()
				ns[i], sk[i] = float(datum[0]), float(datum[1])
			out.close()
		print >> sys.stderr, '\r', "Number of sk configurations:", str(isk), ' ',
		#sys.stderr.flush()

		isk += 1

		#store exponential factors
		for i in xrange(len(x)):
			xbin = int(x[i]*(kmax+1))
			ybin = int(y[i]*(kmax+1))
			zbin = int(z[i]*(kmax+1))
			histo[xbin][ybin][zbin] += 1
#......calculate FFT image of the distribution
		cssf = fft.fftn(histo)

		for l in range(kmax+1):
			for m in range(kmax+1):
				for n in range(kmax+1):
		#tests on magnitude of k vector
					kk = l * l + m * m + n * n

					if kk < ksqmax and kk > 0.:
		#  form exp(ikr) for each particle
		#  form sums for each species
						sk[kk] += abs(cssf[l][m][n]) * abs(cssf[l][m][n])
						ns[kk] += 1
#----------------------------------------------------------------------
	#writing output
 	out = open('sk-' + pid + '.tmp', 'wb')
	s = str(isk) + '\n'
	out.write(s)
	for i in xrange(ksqmax):
 		s =  str(ns[i]) + ' ' + str(sk[i]) + '\n'
		out.write(s)
	out.close()
#--------END OF FUNCTION SK--------------------------------------------

def skfinalize(config, fname):
	pid = str(os.getpid())
	natom = config.n_atom()
	box3d = config.box()
	box = min(box3d.x, box3d.y, box3d.z)
	nmax, kmax, ksqmax = natom, 25, 1024
	ns = zeros(ksqmax, float)
	sk = zeros(ksqmax, float)
	out = open('sk-' + pid + '.tmp', 'r+')
	isk = int(out.readline().split()[0])
 	for i in xrange(ksqmax):
		datum = out.readline().split()
		ns[i], sk[i] = float(datum[0]), float(datum[1])
	out.close()

	outdata = open(fname, 'wb')
	for i in xrange(ksqmax):
		if ns[i] != 0:
			ski = sk[i] / float(natom * natom) / float(ns[i])
			s = str(2.*math.pi*math.sqrt(float(i))/box) + ' ' + str(ski) + '\n'
			outdata.write(s)

	os.remove('sk-' + pid + '.tmp')
	outdata.close()
