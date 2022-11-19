#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 20:15:18 2022

@author: misha
"""

import numpy as np
from neighbor import neighbor_mask

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
        
        out_of_range = neighbor_mask(atom_coordinates, ref_coordinates,
                                      box, r_max = r_max)
        
        masked_neighbors = np.ma.array(atom_indices, mask = out_of_range)
        
        compressed_neighbors = np.ma.compressed(masked_neighbors)
        
        neighborlists[ref_index] = np.ndarray.tolist(compressed_neighbors)
        
    return neighborlists