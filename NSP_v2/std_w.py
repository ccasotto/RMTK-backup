# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 14:11:06 2014

@author: chiaracasotto
"""
import numpy as np
import math

def std_w(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, axis = 0, weights=weights)
    variance = np.average(np.power((values-average),2), axis = 0, weights=weights)  # Fast and numerically precise
    return np.sqrt(variance)