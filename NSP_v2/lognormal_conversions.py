# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 09:54:51 2014

@author: chiaracasotto
"""
import numpy as np

def from_median_to_mean(Sa50,beta):
    # INPUT: Sa50 is the median of x, beta is the disperison of x
    # OUTPUT: Samean is the median of x, StdSa is the std of x
    # Converting medians to means
    Samean = Sa50*np.exp(np.power(beta,2)/2)
    # Converting dispersion std(log(Sa)) into standard deviation std(Sa)
    StdSa = np.sqrt(np.exp(np.power(beta,2)-1)*np.exp(2*np.log(Sa50)+np.power(beta,2)))
    return [Samean,StdSa]

def from_mean_to_median(Samean, StdSa):
    # INPUT: Samean is the median of x, StdSa is the std of x
    # OUTPUT: Sa50 is the median of x, beta is the disperison of x
    Sa50 = np.power(Samean,2)/np.sqrt(np.power(StdSa,2)+np.power(Samean,2)) # median
    beta = np.sqrt(np.log(1+np.power(np.array(StdSa)/np.array(Samean),2))) # dispersion
    return [Sa50, beta]