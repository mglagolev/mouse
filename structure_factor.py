#!/usr/bin/env python3

from lammpack_misc import *
from lammpack_types import *
from ordering_functions import createNumpyAtomsArrayFromConfig
import vector3d
import numpy as np
import sys
import math

def calculateStructureFactor(config, g, ikkmax):
    """ Calculate S values for a range of k^2 values
    Output: S(k^2), number_of_samples(k^2) """
    tbox = 2. * np.pi / float(min(config.box().x, config.box().y, config.box().z))
    sk = np.zeros(ikkmax + 1, dtype = np.longdouble)
    ns = np.zeros(ikkmax + 1, dtype = int)
    npAtoms = createNumpyAtomsArrayFromConfig(config)
    el = np.cos(tbox * npAtoms[3], dtype = np.longdouble) + 1j * np.sin(tbox * npAtoms[3], dtype = np.longdouble)
    em = np.cos(tbox * npAtoms[4], dtype = np.longdouble) + 1j * np.sin(tbox * npAtoms[4], dtype = np.longdouble)
    en = np.cos(tbox * npAtoms[5], dtype = np.longdouble) + 1j * np.sin(tbox * npAtoms[5], dtype = np.longdouble)
    els, ems, ens, ems_neg, ens_neg = [], [], [], [], []
    eln, emn, enn = el[:], em[:], en[:]
    els.append(np.array([1.] * config.n_atom(), dtype = np.clongdouble))
    ems.append(np.array([1.] * config.n_atom(), dtype = np.clongdouble))
    ens.append(np.array([1.] * config.n_atom(), dtype = np.clongdouble))
    for l in range(1, g + 1):
        els.append(eln)
        eln = eln * el
    for m in range(1, g + 1):
        ems.append(emn)
        emn = emn * em
    for n in range(1, g + 1):
        ens.append(enn)
        enn = enn * en
    for m in range(len(ems)):
        ems_neg.append(np.conjugate(ems[m]))
    for n in range(len(ens)):
        ens_neg.append(np.conjugate(ens[n]))
    for l in range(0, g + 1):
        for m in range(-g, g + 1):
            for n in range(-g, g + 1):
                ikk = l * l + m * m + n * n
                if ikk <= ikkmax:
                    if m < 0: em_value = ems_neg[-m]
                    else: em_value = ems[m]
                    if n < 0: en_value = ens_neg[-n]
                    else: en_value = ens[n]
                    expikr = np.sum(els[l] * em_value * en_value)
                    sk[ikk] += np.real(expikr * np.conjugate(expikr))
                    ns[ikk] += 1
    for ikk in range(1, ikkmax+1):
        if ns[ikk] > 0:
            sk[ikk] = sk[ikk] / config.n_atom()
    return [1.] + list(sk[1:]), [1] + list(ns[1:])
