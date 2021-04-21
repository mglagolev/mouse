#!/usr/bin/env python3

from lammpack_misc import *
from lammpack_types import *
from ordering_functions import createNumpyBondsArrayFromConfig
from ordering_functions import createNumpyAtomsArrayFromConfig
import vector3d
import numpy as np
from numpy import linalg as LA
import sys
import pdb

def ConnectingVectorNonnorm(config, com_block_1, com_block_2):
    X_size = config.box().x
    Y_size = config.box().y
    Z_size = config.box().z
    x1 = com_block_1.x
    x2 = com_block_2.x
    y1 = com_block_1.y
    y2 = com_block_2.y
    z1 = com_block_1.z
    z2 = com_block_2.z
    dx = np.remainder(x2 - x1 + float(X_size)/2., float(X_size)) - float(X_size)/2.
    dy = np.remainder(y2 - y1 + float(Y_size)/2., float(Y_size)) - float(Y_size)/2.
    dz = np.remainder(z2 - z1 + float(Z_size)/2., float(Z_size)) - float(Z_size)/2.
    return vector3d.Vector3d(dx, dy, dz)
    
def HistoVectMod(config, diff_array):
    r_mod = []
    h = [0] * 115
    for i in range(0,len(diff_array)):
        dx = diff_array[i].x
        dy = diff_array[i].y
        dz = diff_array[i].z
        r_mod.append(np.sqrt(dx * dx + dy * dy + dz * dz))
    for j in range(0,115):
        d1 = 0.2*j
        d2 = 0.2*(j+1)
        h.append(0)
        for k in range(0,len(diff_array)):
            if r_mod[k]>=d1 and r_mod[k]<d2:
                h[j] += 1
    return h

def EndToEndHist(config, type1 = "1", type2 = "2"):
    # flexible block - atoms belonging to one molecule, type 1
    # helical block - atoms belonging to one molecule, type 2
    used_mol = []
    v_ee_1 = []
    v_ee_2 = []
    v_ee_total = []
    for atom in config.atoms():
      mol_id = atom.mol_id
      if mol_id not in used_mol:
          block_1 = []
          block_2 = []
          for atom in config.atoms():
              if atom.mol_id == mol_id and atom.type == type1:
                  block_1.append(atom)
              if atom.mol_id == mol_id and atom.type == type2:
                  block_2.append(atom)
          used_mol.append(mol_id)
          # find the beginning and end of the blocks
          block_1.sort(key=lambda x: x.num)
          block_2.sort(key=lambda x: x.num)
          start_1 = vector3d.Vector3d(block_1[0].pos.x, block_1[0].pos.y, block_1[0].pos.z)
          fin_1 = vector3d.Vector3d(block_1[-1].pos.x, block_1[-1].pos.y, block_1[-1].pos.z)
          start_2 = vector3d.Vector3d(block_2[0].pos.x, block_2[0].pos.y, block_2[0].pos.z)
          fin_2 = vector3d.Vector3d(block_2[-1].pos.x, block_2[-1].pos.y, block_2[-1].pos.z)
          v_ee_1.append(ConnectingVectorNonnorm(config, start_1, fin_1))
          v_ee_2.append(ConnectingVectorNonnorm(config, start_2, fin_2))
          v_ee_total.append(ConnectingVectorNonnorm(config, start_1, fin_2))
     
    h_1 = HistoVectMod(config, v_ee_1)
    h_2 = HistoVectMod(config, v_ee_2)
    h_total = HistoVectMod(config, v_ee_total)
    
    for j in range(0,115):
         sys.stderr.write(str(0.2*j+0.1)+" "+str(h_1[j])+" "+str(h_2[j])+" "+str(h_total[j])+"\n")
          
