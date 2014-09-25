# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 12:26:31 2014

@author: chiaracasotto
"""

import numpy as np
import scipy.stats as stat
import os

def damage_to_loss(SaT,bTSa,iml):
    cd = os.getcwd()
    # INPUT: SaT is the mean log(iml), bTSa is the dispersion, std(log(iml))
    # Input consequence functions. The number of damage state considered should
    # be the same in both fragility and consequence. 
    input2 = cd+'/inputs/consequence.csv'
    with open(input2, 'rb') as f:
        data = f.read()
        l = data.rstrip()
        lines = l.split('\n')   
        newlist = [lines[i].split(',') for i in range(0, len(lines))]   
        loss_ratio=[float(ele) for ele in newlist[1]]

    # From continuous function to dicrete
    poe = [np.zeros_like(iml)]*len(SaT)
    for i in range(0,len(SaT)):
        mu = SaT[i]
        sigma = bTSa[i]
        if sigma == 0:
            probability = np.zeros_like(iml)
            probability[np.log(iml)>mu] = 1 # check this
            poe[i] = probability
        else:
            poe[i] = stat.norm(mu,sigma).cdf(np.log(iml)) # Probability of exceedance each DS
    
    # From probability of exceedance to probability of occurance and loss ratio per damage state
    loss_ratio_DS = [np.zeros_like(iml)]*len(SaT)
    poo = poe
    for i in range(0,len(SaT)-1):
        poo[i] = poe[i]-poe[i+1] # probability of occurance of each DS
    for i in range(0,len(SaT)):    
        loss_ratio_DS[i] = poo[i]*loss_ratio[i] # loss_ratio vs intensity measure for each DS

    # Total loss_ratio
    LR = sum(loss_ratio_DS)
    return LR