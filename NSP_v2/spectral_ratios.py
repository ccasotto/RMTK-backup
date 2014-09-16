# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 14:40:47 2014

@author: chiaracasotto
"""
import os
import numpy as np

def get_spectral_ratios(Tuni,T):
    cd = os.getcwd()
    input = cd+'/inputs/FEMAP965spectrum.txt'
    with open(input, 'rb') as f:
            data = f.read()
            l = data.rstrip()
            lines = l.split('\n')
            data = [lines[i].split('\t') for i in range(0, len(lines))]
            Ts = np.array([float(ele[0]) for ele in data])
            Sa = np.array([float(ele[1]) for ele in data])
            S_Tuni = np.interp(Tuni,Ts,Sa)
            S_T = np.array([np.interp(ele,Ts,Sa) for ele in T])
            Sa_ratios = S_Tuni/S_T
            
    return Sa_ratios