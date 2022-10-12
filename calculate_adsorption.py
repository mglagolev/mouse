#!/usr/bin/env python3

import sys
import os
import argparse
import math
import statistics
import string
import vector3d
import lammpack_misc
from lammpack_types import Config
import numpy as np

parser = argparse.ArgumentParser(description = 'Calculate msd of the adsorbed molecules')

parser.add_argument('datafiles', metavar = 'LAMMPS_data', type = str, nargs = '+', help = 'lammps data file')

parser.add_argument('--dumps', type = str, nargs = '+', help = 'additional dump files')

parser.add_argument('--nchains', type = int, nargs = '?', help = 'number of chains')

parser.add_argument('--nlen', type = int, nargs = '?', help = 'chain length')

parser.add_argument('--ntotal', type = int, nargs = '?', help = 'total number of atoms')

parser.add_argument('--zsorb', type = float, nargs = '?', help = 'z for adsorbed atoms')

args = parser.parse_args()

#print(args.zsorb)

if len(args.datafiles) > 1 and len(args.dumps) > 0:
    raise NameError("Additional dump files can be used only with a single data file")

configs = []
for in_data in args.datafiles:
    config = Config()
    config.read_lmp_data(in_data)
    configs = [config]
    if args.dumps is not None:
        for dumpfile in args.dumps:
            extra_configs = config.return_extra_frames_from_lmp_dump(dumpfile)
            extra_configs_exclude_datafile = [x for x in extra_configs if x.timestep != config.timestep]
            configs += extra_configs_exclude_datafile
            
timesteps = [x.timestep for x in configs]
timesteps.sort()

print(timesteps)

def fix_uncut_coords(configs):
    """ Consider coordinates in the first config as uncut,
    for each of the next configs calculate the uncut coordinates
    as the image closest to the previous coordinates"""
    atom_ids = [x.id for x in configs[0].atoms()]
    for atom_id in atom_ids:
        for i in range(len(configs)):
            config = configs[i]
            atom = config.atom_by_id(atom_id)
            if i == 0:
                x_uncut_cur = atom.pos.x + config.box().x * atom.pbc.x
                y_uncut_cur = atom.pos.y + config.box().y * atom.pbc.y
                z_uncut_cur = atom.pos.z + config.box().z * atom.pbc.z
                pos_uncut_cur = vector3d.Vector3d(x_uncut_cur, y_uncut_cur, z_uncut_cur)
            else:
                shift = atom.pos - pos_uncut_cur #pos_cut_new - pos_cur
                shift_trimmed = lammpack_misc.vector_pbc_trim(shift, config.box())
                pos_uncut_new = pos_uncut_cur + shift_trimmed
                ix = math.floor((pos_uncut_new.x - config.box_center().x + config.box().x / 2. ) / config.box().x)
                iy = math.floor((pos_uncut_new.y - config.box_center().y + config.box().y / 2. ) / config.box().y)
                iz = math.floor((pos_uncut_new.z - config.box_center().z + config.box().z / 2. ) / config.box().z)
                atom.pbc = vector3d.Vector3d(ix, iy, iz)
                pos_uncut_cur = pos_uncut_new            
            
fix_uncut_coords(configs)

adsorbtion_data = {}


mol_ids = [x.mol_id for x in config.atoms() if x.mol_id != "0" ]
unique_mol_ids = set(mol_ids)

for mol_id in unique_mol_ids:
    for timestep in timesteps:
        config = [x for x in configs if x.timestep == timestep][0]
        molecule_atoms = [x for x in config.atoms() if x.mol_id == mol_id]
        is_adsorbed = False
        for atom in molecule_atoms:
            if atom.pos.z <= args.zsorb and (atom.type == '3' or atom.type == '4'):
                is_adsorbed = True
        try:
            adsorbtion_data[mol_id].append([timestep, is_adsorbed])
        except KeyError:
            adsorbtion_data[mol_id] = [[timestep,is_adsorbed]]

            
tracks = []
for mol_id in unique_mol_ids:
    """ Write down each interval when the molecule is adsorbed into
    an individual track, and add the molecule id. The format is:
    [mol_id, [timestep0, timestep1, ...]]"""
    is_molecule_adsorbed = adsorbtion_data[mol_id]
    state = False
    for i in range(len(is_molecule_adsorbed)):
        currently_adsorbed = is_molecule_adsorbed[i][1]
        timestep = is_molecule_adsorbed[i][0]
        if state == False and currently_adsorbed == True:
            track = [timestep]
            state = True
        elif state == True and currently_adsorbed == True:
            track.append(timestep)
            if len(track) == len(timesteps):
                tracks.append([mol_id,track])
        elif state == True and currently_adsorbed == False:
            if len(track) >= 2:
                tracks.append([mol_id,track])
                state = False
        elif state == False and currently_adsorbed == False:
            pass

def calculate_msd_for_track(mol_id, track, configs):
    """ Calculate the mean-square displacement of the center of mass
    relative to the first timestep of the track"""
    msd = 0.0
    t = 0
    print(mol_id)
    msd_track = []
    for i in range(len(track[1])):
        step = track[1][i]
        print('step',step)
        config = [x for x in configs if x.timestep == step][0]
        molecule_atoms = [x for x in config.atoms() if x.mol_id == mol_id]
        sum_X = 0.0
        sum_Y = 0.0
        sum_Z = 0.0 
        for atom in molecule_atoms:
            sum_X += config.box_center().x + atom.pos.x + config.box().x * (atom.pbc.x - 0.5)
            sum_Y += config.box_center().y + atom.pos.y + config.box().y * (atom.pbc.y - 0.5)
            sum_Z += config.box_center().z + atom.pos.z + config.box().z * (atom.pbc.z - 0.5)
        X_com1 = sum_X/float(len(molecule_atoms))
        Y_com1 = sum_Y/float(len(molecule_atoms))
        Z_com1 = sum_Z/float(len(molecule_atoms))
        if t == 0:
            x0 = X_com1
            y0 = Y_com1
            z0 = Z_com1
        if t != 0:
            msd = (X_com1 - x0) * (X_com1 - x0) + (Y_com1 - y0) * (Y_com1 - y0)
        msd_track.append(msd)
        t = t + 1
    return msd_track    

msds_data = []

M = 0    
f_2 = open("out-tracks.txt","w+")
print(tracks)                
for track in tracks:
    """Calculate MSDs from each track"""
    print(track)
    mol_id = track[0]
    timesteps = track[1]
    for ele in track:
        f_2.write(str(ele)+"\n")
    print(track)
    msds = calculate_msd_for_track(mol_id, track, configs)
    if len(msds)>M:
    	M = len(msds)
    for ele in msds:
        f_2.write(str(mol_id)+" "+str(ele)+"\n")
    print('mol_id=', mol_id, 'msds=', msds)
    try:
        msds_data.append(msds)
    except KeyError:
        msds_data = msds
f_2.close()
    
f_1 = open("out-ave-msd.txt","w+")
max_msd_len = max(map(len,msds_data))
ave_msd = []
for i in range(max_msd_len):
    data = [ x[i] for x in msds_data if len(x) > i ]
    ave_msd.append(sum(map(float,data))/len(data))
    f_1.write("%d %7.3f \n" % (i, ave_msd[i]))
    print(i, ave_msd[i])
f_1.close()
