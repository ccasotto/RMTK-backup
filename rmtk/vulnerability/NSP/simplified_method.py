# -*- coding: utf-8 -*-
"""
Created on Thu May 29 11:29:32 2014

@author: chiaracasotto
"""
import numpy as np
import scipy.stats as stat
import os
pi = 3.141592653589793

def simplified_bilinear(T, Gamma, dcroof, SPO, bUthd, g):
    # Step 1: Ductility level mu for each Limit State
    dy = SPO[0]
    mu = np.divide(dcroof,dy)
    print "mu(LS) = ", mu
    # Step 2: Define constant relative strength inelastic displacement ratio Cr
    c = np.multiply(79.12,np.power(T,1.98))
    R1 = 0.425*(1-c+np.sqrt(np.power(c,2)+c*2*(2*mu-1)+1))
    R = np.array([max(ele,1.00001) for ele in R1])
    Cr50 = 1+np.divide(R-1,c)

    # Step 4: Estimate Sa50 (median IM per Limit State) in m/s^2
    Sa50 = np.divide(np.power(2*pi/T,2),g*Gamma*Cr50)*dcroof 
    # Step 5: Estimate b
    b = 1+np.divide(np.log(Cr50),np.log(R))

    # Step 6: Estimate dispersion of theta demand
    sigmalnd = 1.957*(1/5.876+1/(11.749*(T+0.1)))*(1-np.exp(-0.73*(R-1)))
    bthd = sigmalnd
    bTSa = 1/b*np.sqrt(np.power(bthd,2)+np.power(bUthd,2))
    print "median IM = ", Sa50
    print "total dispersion = ", bTSa
    return [Sa50, bTSa]