#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 17:20:41 2022

@author: Mikhail Glagolev
"""

import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np
from ordering import calculate_cos_sq_for_reference

def calculate_orientation_order_parameter(
        u: mda.Universe, r_min = 0., r_max = -1., n_bins = 150,
        mode = 'average', same_molecule = True):
    """
    
    This function calculates local ordering parameter for bonds.
    The parameter is calculated as S = 3/2 ( <(cos(gamma))^2>) - 1/2)
    where "gamma" is the angle between the bond vectors. The distance
    between the middle points of the bond vectors shall be in the range
    [rmin, rmax].

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis universe. Only the bonds required for calculation of the
        ordering parameter shall be present. All other bonds shall be deleted
        from the universe before the analysis.
    r_min : FLOAT, optional
        Minimum distance between the bond vector centers to consider.
        The default is 0.. To exclude the bond itself, consider setting
        r_min to a small value, e. g. 1e-6
    r_max : FLOAT, optional
        Maximum distance between the bond vector centers to consider.
        The default is 0., which means choosing the cutoff based on the
        size of the simulation cell.
        Setting the value to -1, means considering all the bonds.
    mode : STRING, optional
        Whether an average value or a histogram shall be returned.
        The default is 'average'.
    same_molecule : BOOL, optional
        Whether the bonds from the same molecule (resid) shall be accounted
        for. The default is True.

    Returns
    -------
    FLOAT
        Average value of the local ordering parameter.

    """
    if mode == "average":
        cos_sq_sum = 0.
        i_s = 0
    elif mode == "histogram":
        cos_sq_raw_hist = np.zeros(n_bins)
        _, bin_edges = np.histogram(cos_sq_raw_hist, bins = n_bins,
                                           range = (0.,1.))
        
    if r_max == 0.:
        r_max = min(u.dimensions) / 2.
    # Unwrap all the coordinates, so that all the bond lengths are
    # real. The closest images of the bonds will be found in the nested
    # function.
    unwrap = transformations.unwrap(u.atoms)
    u.trajectory.add_transformations(unwrap)
    
    # Calculate bond components
    # 1D arrays, one for each of the coordinates, provide more efficient
    # numpy calculations. Converting the data here, outside of the main loop
    # provided additional 15% speedup in the test runs.
    bx = (u.bonds.atom2.positions[:, 0] - u.bonds.atom1.positions[:, 0])
    by = (u.bonds.atom2.positions[:, 1] - u.bonds.atom1.positions[:, 1])
    bz = (u.bonds.atom2.positions[:, 2] - u.bonds.atom1.positions[:, 2])
    
    bond_components = [bx, by, bz]
    
    # Creating 1D arrays with bond midpoints
    rx = (u.bonds.atom1.positions[:, 0] + u.bonds.atom2.positions[:, 0]) / 2.
    ry = (u.bonds.atom1.positions[:, 1] + u.bonds.atom2.positions[:, 1]) / 2.
    rz = (u.bonds.atom1.positions[:, 2] + u.bonds.atom2.positions[:, 2]) / 2.
    
    bond_midpoints = [rx, ry, rz]
    
    if not same_molecule:
        bond_resids = u.bonds.atom1.resids
    else:
        bond_resids = None
    
    for bond in u.bonds:
        # Determine the reference vector components and midpoint
        # from the bond coordinates
        ref_components = bond.atoms[1].position - bond.atoms[0].position
        ref_midpoint = (bond.atoms[0].position + bond.atoms[1].position) / 2.
        # If needed, exclude bonds from the same molecule
        if not same_molecule:
            excluded_resids = bond.atoms[0].resid
        else:
            excluded_resids = None
        # Calculate ordering parameter value for the reference bond
        cos_sq_masked = calculate_cos_sq_for_reference(
            bond_components, bond_midpoints, ref_components, ref_midpoint,
            u.dimensions, r_min = r_min, r_max = r_max,
            vector_attributes = bond_resids,
            excluded_attributes = excluded_resids)
            
        if mode == "average":
            if np.ma.count(cos_sq_masked) > 0:
                cos_sq_sum += np.ma.average(cos_sq_masked)
                i_s += 1
            else:
                pass
                
        elif mode == "histogram":
            cos_sq_hist_increment, _ = np.histogram(
                np.ma.compressed(cos_sq_masked),
                bins = n_bins, range = (0.,1.))
            cos_sq_raw_hist += cos_sq_hist_increment

    if mode == "average":
        if i_s > 0:
            # Normalize the values. Normalization procedure ensures that double
            # consideration of each of the bonds doesn't affect the result
            s = 1.5 * cos_sq_sum / i_s - 0.5
            return s
        else:
            raise NameError("No pairs of bonds found within range")
    elif mode == "histogram":
        norm = np.sum(cos_sq_raw_hist * np.diff(bin_edges))
        print("Norm " + str(norm) + "\n")
        cos_sq_hist = cos_sq_raw_hist / norm
        return cos_sq_hist

if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser(
        description = 'Calculate local ordering parameter')

    parser.add_argument(
        'input', metavar = 'INPUT', action = "store", help = "input file")

    parser.add_argument(
        '--r_max', metavar = 'R_max', type = float, nargs = '?',
        default = 0., help = "outer cutoff radius")

    parser.add_argument(
        '--r_min', metavar = 'R_min', type = float, nargs = '?',
        default = 0., help = "inner cutoff radius")

    parser.add_argument(
        '--mode', metavar = 'MODE', type = str, nargs = '?',
        default = "average", help = "mode: average or historgram")

    parser.add_argument(
        '--n_bins', metavar = 'N_bins', type = int, nargs = '?',
        default = 150, help = "number of bins of the histogram")

    parser.add_argument(
        "--same-molecule", action = "store_true",
        help = "Hide the bonds transversing the periodic boundary conditions")

    args = parser.parse_args()

    u = mda.Universe(args.input)
    
    s = calculate_orientation_order_parameter(u, r_min = args.r_min,
                                              r_max = args.r_max,
                                              mode = args.mode,
                                              n_bins = args.n_bins,
                                              same_molecule
                                              = args.same_molecule)

    print(str(s) + "\n")