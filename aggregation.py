#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 20:15:18 2022

@author: misha
"""

import numpy as np

def calculate_neighborlists_from_distances(
        atom_indices: np.ndarray,
        atom_coordinates: [np.ndarray, np.ndarray, np.ndarray],
        box = [0., 0., 0.], r_max = 0.):
    
    neighborlists = {}
    
    if (min(atom_indices.size, atom_coordinates[0].size,
           atom_coordinates[1].size, atom_coordinates[2].size) !=
       max(atom_indices.size, atom_coordinates[0].size,
              atom_coordinates[1].size, atom_coordinates[2].size)):
        raise NameError("The dimensions of input datasets differ")
    
    #Iterate over atoms
    for i in range(atom_indices.size):
        
        ref_index = atom_indices[i]
        
        ref_coordinates = [atom_coordinates[0][i],
                           atom_coordinates[1][i],
                           atom_coordinates[2][i]]
        
        #TODO: this section is identical to the corresponding one in ordering.py
        #Consider creating a specialized function for finding the closest neighbor
        #TODO: And check the performance
        
        drx = atom_coordinates[0] - ref_coordinates[0]
        dry = atom_coordinates[1] - ref_coordinates[1]
        drz = atom_coordinates[2] - ref_coordinates[2]

        #Find the bond image which is closest to the reference vector
        if box[0] > 0.:
            drx = drx - box[0] * np.around(drx / box[0])
        if box[1] > 0.:
            dry = dry - box[1] * np.around(dry / box[1])
        if box[2] > 0.:
            drz = drz - box[2] * np.around(drz / box[2])

        #Calculate the distance between the atoms
        rsq = drx**2 + dry**2 + drz**2
        
        out_of_range = np.greater(rsq, r_max**2)
        
        masked_neighbors = np.ma.array(atom_indices, mask = out_of_range)
        
        compressed_neighbors = np.ma.compressed(masked_neighbors)
        
        neighborlists[ref_index] = np.ndarray.tolist(compressed_neighbors)
        
    return neighborlists