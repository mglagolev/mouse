#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 14:46:18 2022

@author: Mikhail Glagolev
"""

import numpy as np
from neighbor import neighbor_mask

    
def calculate_cos_sq_for_reference(
        vector_components: [np.ndarray, np.ndarray, np.ndarray],
        vector_midpoints: [np.ndarray, np.ndarray, np.ndarray],
        ref_components: np.ndarray(3), ref_midpoint: np.ndarray(3),
        box = [0., 0., 0., 90., 90., 90.], r_min = 0., r_max = -1.,
        vector_attributes = None, excluded_attributes = None):
    """
    

    Parameters
    ----------
    vector_components : [np.ndarray, np.ndarray, np.ndarray]
        Components of the vectors in Cartesian coordinates
    vector_midpoints : [np.ndarray, np.ndarray, np.ndarray]
        Midpoints of the vectors in Cartesian coordinates
    ref_components : np.ndarray(3)
        Components of the reference vector in Cartesian coordinates
    ref_midpoint : np.ndarray(3)
        Midpoint of the reference vector in Cartesian coordinates
    box : TYPE, optional
        If the periodic boundary conditions are to be taken into account,
        the box dimensions shall be provided. Along the axes where the
        dimensions of the box are 0, the periodic boundary conditions are
        not taken into account. The default is [0., 0., 0.].
    r_min : TYPE, optional
        Internal cutoff radius for the distance between the vector midpoints.
        The default is 0..
    r_max : TYPE, optional
        External cutoff radius for the distance between the vector midpoints.
        If 0. is provided, the value will be determined as half of the lowest
        of the dimensions of a periodic box.
        The default is -1, corresponding to no cutoff
    vector_attributes : np.ndarray, optional
        An array of attributes associated with the vectors. This can be used
        if some vectors shall not be accounted for. The default is None.
    excluded_attributes : optional
        The values of the attribute, when the vectors shall not be accounted
        for. The default is None.

    Raises
    ------
    NameError
        "No values to calculate"
        If no pairs of vectors were found within the cutoff distance, this
        error will be raised.

    Returns
    -------
    numpy.ma.masked_array
        The array with the resulting values together with the mask, marking
        the values that do not satisfy the selection criteria: cutoff radii
        and blacklisted attributes

    """
    # If r_max is set to 0, determine the largest possible cutoff based on box.
    if r_max == 0.:
        r_max = np.min(box / 2.)
    
    if np.linalg.norm(ref_components) == 0.:
        raise NameError("Zero length of reference vector")

    # The backend can be "NumPy", "MDA-serial", "MDA-OpenMP"
    out_of_range = neighbor_mask(vector_midpoints, ref_midpoint,
                                      box, r_min, r_max, backend = "NumPy")
    
    # Create an array masking the values that shall be excluded
    if excluded_attributes is not None:
        excluded = np.equal(vector_attributes, excluded_attributes)
    else:
        excluded = np.zeros(vector_components[0].size)
    
    masked_data = np.logical_or(out_of_range, excluded)
    
    # Calculate the nparray for normalized cos^2(theta), where theta are the
    # angles between the bonds and the reference vector
    cos_sq_normed = (
            (ref_components[0] * vector_components[0]
           + ref_components[1] * vector_components[1]
           + ref_components[2] * vector_components[2] )**2
          / np.linalg.norm(ref_components)**2
          / np.linalg.norm(vector_components, axis = 0)**2)
    
    # Mask the invalid values and return
    return np.ma.array(cos_sq_normed, mask = masked_data)