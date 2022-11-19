#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 23:34:30 2022

@author: Mikhail Glagolev

TODOs:

 - support working with multiple timesteps 
 (mind that the box size can vary between the timesteps)
"""

import MDAnalysis as mda
import numpy as np
import networkx as nx
from aggregation import calculate_neighborlists_from_distances

def determine_aggregates(u: mda.Universe, r_max: float, selection = None):
    # Select atoms by type or read selection criteria in MDAnalysis synthax
    if selection is not None:
        atoms = u.select_atoms(selection)
    else:
        atoms = u.atoms
    #A list of lists [[aggregate1_atom1, aggregate1_atom2, ...],
    # [aggregate2_atom1, aggregate2_atom2, ...], [aggregate3_atom1, ...]]
    aggregates = []
    # Create numpy array with a list of atom indices
    atom_indices = atoms.indices
    # Create numpy array with the coordinates of the atoms
    rx = atoms.positions[:, 0]
    ry = atoms.positions[:, 1]
    rz = atoms.positions[:, 2]
    
    atom_positions = [rx, ry, rz]
    # Create neighbor lists. Compare my function to the MDAnalysis standard
    # function
    # neighborlists = { atom_index : [neighbor1_index, neighbor2_index, ...]}
    neighborlists = calculate_neighborlists_from_distances(atom_indices,
                                                           atom_positions,
                                                           box = u.dimensions,
                                                           r_max = r_max)
    #Option 2: use MDAnalysis function to calculate neighbor lists:
    
    # Initialize a NetworkX graph
    graph = nx.Graph()
    # For every atom add the neighbors to the graph
    for atom_index in neighborlists:
        for neighbor in neighborlists[atom_index]:
            graph.add_edge(atom_index, neighbor)
    # Convert atom indices to a list
    atom_indices_list = np.ndarray.tolist(atom_indices)
    # While list length > 0:
    while len(atom_indices_list) > 0:
        aggregate = []
        aggregate_atoms = nx.dfs_postorder_nodes(graph, atom_indices_list[0])
        for atom in aggregate_atoms:
            aggregate.append(atom)
            atom_indices_list.remove(atom)
        aggregates.append(aggregate)
    return aggregates

if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser(
        description = 'Determine aggregates, based on lists of neighbors\
        determined by inter-particle distance')

    parser.add_argument(
        'input', metavar = 'INPUT', action = "store", help = "input file")

    parser.add_argument(
        '--r_max', metavar = 'R_max', type = float, nargs = '?',
        default = 1.2, help = "neighbor cutoff")

    parser.add_argument(
        '--selection', metavar = 'QUERY', type = str, nargs = '?',
        help 
        = "Consider only selected atoms, use MDAnalysis selection language")

    args = parser.parse_args()

    u = mda.Universe(args.input)

    determine_aggregates(u, args.r_max, args.selection)