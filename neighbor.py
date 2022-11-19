#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 13:36:56 2022

@author: Mikhail Glagolev
"""

import numpy as np

def calculate_squared_distances(coordinates, reference_coordinates, box):
    """
    This function employs numpy to calculate an array of squared distances
    between the points with the coordinates defined by three input arrays,
    [x_coordinates, y_coordinates, z_coordinates] and the reference point
    with the coordinates defined by reference_coordinates=[x_ref, y_ref, z_ref]
    
    The periodic boundary conditions will be taken into account for the
    dimensions where the box size is non-zero.
    
    The function returns a numpy array of squared distances.

    Parameters
    ----------
    coordinates : [np.ndarray, np.ndarray, np.ndarray]
        DESCRIPTION.
        
    reference_coordinates : np.ndarray(3)
        The cartesian coordinates of a reference point.
        
    box : np.ndarray, with shape[0] >= 3.
        The box from MDAnalysis library,
        which contains dimensions and angles, can be provided as is
        

    Returns
    -------
    rsq: np.ndarray
        A numpy array of squared distances between the points with coordinates
        defined by the input arrays and the reference point.

    """
    drx = coordinates[0] - reference_coordinates[0]
    dry = coordinates[1] - reference_coordinates[1]
    drz = coordinates[2] - reference_coordinates[2]

    #Find the bond image which is closest to the reference vector
    if box[0] > 0.:
        drx = drx - box[0] * np.around(drx / box[0])
    if box[1] > 0.:
        dry = dry - box[1] * np.around(dry / box[1])
    if box[2] > 0.:
        drz = drz - box[2] * np.around(drz / box[2])
        
    #Calculate the distance between the atoms
    rsq = drx**2 + dry**2 + drz**2
    
    return rsq

def neighbor_mask(coordinates: [np.ndarray, np.ndarray, np.ndarray],
                  reference_coordinates: np.ndarray(3),
                  box, r_min = 0., r_max = 0.):
    """
    
    This function returns the mask which can be applied to a numpy array
    to mask the points that are not neighbors of the reference point.
    The neighborhood is determined by criterion r_min <= r <= r_max, where
    r is the distance between the point determined by the coordinates in the 
    input arrays and the reference point. For non-zero components of box,
    the periodic boundary conditions are taken into account.

    Parameters
    ----------
    coordinates : [np.ndarray, np.ndarray, np.ndarray]
        DESCRIPTION.
        
    reference_coordinates : np.ndarray(3)
        The cartesian coordinates of a reference point.
        
    box : np.ndarray, with shape[0] >= 3.
        The box from MDAnalysis library,
        which contains dimensions and angles, can be provided as is
        
    r_min: float
        The minimum distance between the points to be considered neighbors
        
    r_max : float
        The maximum distance between the points to be considered neighbors

    Returns
    -------
    out_of_range = np.ndarray
        This is a mask that can be applied to a numpy array to mask the points
        that are not neighbors of the reference point.

    """
    
    rsq = calculate_squared_distances(coordinates, reference_coordinates, box)
    
    # consider r_min as non-negative. For r_min <= 0, the check can be omitted.
    if r_min > 0.:
        # to omit the r_max check, a negative value can be provided, r_max = -1
        if r_max >= 0.:
            out_of_range = np.logical_or(np.less(rsq, r_min**2),
                              np.greater(rsq, r_max**2))
        else:
            out_of_range = np.less(rsq, r_min**2)
    else:
        if r_max >= 0.:
            out_of_range = np.greater(rsq, r_max**2)
        else:
            out_of_range = np.zeros(coordinates[0].size)
    
    return out_of_range