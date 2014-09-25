# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 10:03:12 2014

@author: chiaracasotto
"""
import numpy as np

def count_to_poe(dcm,totblg):
    dpm = np.matrix(np.empty([dcm.shape[0],dcm.shape[1]]))
    fr = np.matrix(np.empty([dcm.shape[0],dcm.shape[1]])) #change here dimensions
    for line in range(0,dcm.shape[0]):
        for ls in range(0,dcm.shape[1]):
            dpm[line,ls] = np.divide(dcm[line,ls],totblg) # probability of occurance of damage
            acc = -dpm[line,:ls]
            cc = acc.tolist()
            cdm = [1]+cc[0]
            fr[line,ls] = sum(cdm) # probability of exceedance of damage
            fr = fr.round(5)
            fr = fr.__abs__()
    return [fr]