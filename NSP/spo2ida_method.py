# -*- coding: utf-8 -*-
"""
Created on Thu May 29 11:29:32 2014

@author: chiaracasotto
"""
import numpy as np
import scipy.stats as stat
import os
pi = 3.141592653589793

def spo2ida(idacm, idacr, mf, T, Gamma, g, dcroof, SPO, bUthd, MC):
    dry = SPO[0]
    mcroof = np.array(dcroof)/(dry) # ductility levels at limit states
    print "mu(LS) = ", mcroof
    
    # Assume lognormal and do some Monte Carlo
    musample = []
    for i in range(0,len(mcroof)):
        musample.append([])
    
    st = (1./(2.*MC))
    en = (1.-(1./(2.*MC)))
    xp = np.linspace(st,en,MC)
    for i in range(0,len(mcroof)):
        if bUthd[i]>0:            
            musample[i] = stat.lognorm.ppf(xp,bUthd[i],loc=0,scale=mcroof[i])
            musample[i][musample[i]>mf]=mf
        else:
            musample[i] = np.repeat(mcroof[i],MC)
    
    Rcap = [[],[],[]]
    RcapMC = [[],[],[]]
    for j in range(0,3):
        for i in range(0,len(mcroof)):
            RcapMC[j].append([])

    # find where in ida curves ductility (i) is reached and get corresponding R for all the percentile (j)
    for i in range(0,len(mcroof)):
        for j in range(0,3):
            Rcap[j].append(np.interp(mcroof[i],idacm[j],idacr[j]))
            RcapMC[j][i]=(np.interp(musample[i],idacm[j],idacr[j]))
                # RcapMC[j][0] with j =0,1,2 are the 3 percentiles corresponding to musample ductilities
                
    Say = np.power(2*pi,2)*dry/(g*Gamma*np.power(T,2))
    # Step 4: Estimate Sa50 (median IM per Limit State) in m/s^2 without considering bUthd
    SaR50 = np.array(Rcap[1])*Say
    bRSa = 0.5*(np.log(np.array(Rcap[0]))-np.log(np.array(Rcap[2])))

#    # These are simple estimates of relatively good accuracy
#    SaT50_0 = [np.median(RcapMC[1][i])*Say for i in range(0,len(mcroof))]
#    bRc = [np.std(np.log(RcapMC[1][i])) for i in range(0,len(mcroof))]
#    bTSa_0 = np.sqrt(np.power(bRSa,2) + np.power(np.array(bRc),2))
    
    Sai = []
    for i in range(0,len(mcroof)):
        Sai.append([])

    # this is to compute the better estimate, still it is much different from above
    # estimates (at least for most cases that I have tried).
    if np.array(bUthd).any()>0:
        allSa50 = [ele*Say for ele in RcapMC[1]]
        allbSa50 = [((np.log(RcapMC[0][h])-np.log(RcapMC[2][h]))/2) for h in range(0,len(mcroof))] #check this        
        for i in range(0,len(mcroof)):      
            for j in range(0,MC):
                if allbSa50[i][j]== 0:
                    realisation = np.repeat(allSa50[i][j],MC)
                    Sai[i].append(realisation)
                else:
                    realisation = stat.lognorm.ppf(xp,allbSa50[i][j],loc=0,scale=allSa50[i][j])
                    Sai[i].append(realisation)                
        SaT50=[]
        bTSa = []
        
        Sa = []
        for i in range(0,len(mcroof)):
            Sa.append(np.array([]))
            
        for i in range(0,len(mcroof)):
            for j in range(1,len(Sai[i])):
                Sai[i][j] = np.concatenate((Sai[i][j-1],Sai[i][j]))
            Sa[i] = Sai[i][-1]
            
        SaT50=[np.median(ele) for ele in Sa]
        bTSa=[np.std(np.log(ele)) for ele in Sa]
    else:
        SaT50 = SaR50
        bTSa = bRSa
    print "median IM = ", SaT50
    print "total dispersion = ", bTSa
    return [SaT50,bTSa]
