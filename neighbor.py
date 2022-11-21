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
        Three 1D arrays, containing the values of the corresponding
        coordinate for all of the points.
        
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

def calculate_distances(coordinates, reference_coordinates, box):
    """
    
    This function calculates distances between the points defined by
    the coordinates array [x_array, y_array, z_array] and the reference
    point [x_ref, y_ref, z_ref], taking into account the periodic boundary
    conditions defined by the box.
    It takes the square root from the squared distances calculated by the
    corresponding function.
    
    Parameters
    ----------
    coordinates : [np.ndarray, np.ndarray, np.ndarray]
        Three 1D arrays, containing the values of the corresponding
        coordinate for all of the points.
        
    reference_coordinates : np.ndarray(3)
        The cartesian coordinates of a reference point.
        
    box : np.ndarray, with shape[0] >= 3.
        The box from MDAnalysis library,
        which contains dimensions and angles, can be provided as is
        

    Returns
    -------
    r: np.ndarray
        A numpy array of distances between the points with coordinates
        defined by the input arrays and the reference point.
    """
    rsq = calculate_squared_distances(coordinates, reference_coordinates, box)
    
    r = np.sqrt(rsq)
    
    return r

def neighbor_mask(coordinates: [np.ndarray, np.ndarray, np.ndarray],
                  reference_coordinates: np.ndarray(3),
                  box, r_min = 0., r_max = 0., backend = "NumPy"):
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
        
    backend: string
        Backend to calculate the distances between the points in the array
        and the reference point.
        Currently supported options are:
            NumPy: use built-in code based on NumPy vector operations
                (default). According to the tests, this is currently the
                fastest option.
            MDA-serial: use the "distance_array" function from MDAnalysis
                with the "serial" backend.
            MDA-OpenMP: use the "distance_array" function from MDAnalysis
                with the "OpenMP" backend.

    Returns
    -------
    out_of_range = np.ndarray
        This is a mask that can be applied to a numpy array to mask the points
        that are not neighbors of the reference point.

    """
    # MDAnalysis currently supports two of it's own backends to calculate the
    # neighbor matrix, "serial" and "OpenMP". This function will take the
    # values "MDA-serial" or "MDA-OpenMP" and use MDAnalysis with the
    # corresponding backend
    if backend[:3] == "MDA":
        import MDAnalysis as mda
        stacked_coordinates = np.column_stack((
            coordinates[0], coordinates[1], coordinates[2]))
        r = mda.lib.distances.distance_array(
            stacked_coordinates, reference_coordinates, box,
            backend = backend[4:])
        values = r.reshape((-1,))
        min_value = r_min
        max_value = r_max
    # Use built-in NumPy function
    elif backend == "NumPy":
        values = calculate_squared_distances(
            coordinates, reference_coordinates, box)
        min_value = r_min**2
        max_value = r_max**2
    
    # consider r_min as non-negative. For r_min <= 0, the check can be omitted.
    if r_min > 0.:
        # to omit the r_max check, a negative value can be provided, r_max = -1
        if r_max >= 0.:
            out_of_range = np.logical_or(np.less(values, min_value),
                              np.greater(values, max_value))
        else:
            out_of_range = np.less(values, min_value)
    else:
        if r_max >= 0.:
            out_of_range = np.greater(values, max_value)
        else:
            out_of_range = np.zeros(coordinates[0].size)
    
    return out_of_range