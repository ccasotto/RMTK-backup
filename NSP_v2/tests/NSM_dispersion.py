# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

"""
Created on Tue May 13 18:04:57 2014

@author: chiaracasotto
"""
# Clear existing variables
def clearall():
    all = [var for var in globals() if var[0] != "_"]
    for var in all:
        del globals()[var]
clearall()

# Import functions
import matplotlib.pyplot as plt
import numpy as np
import os
from common.get_data import read_data
from common.get_Sc_bt import spo2ida
from common.get_Sc_bt import simplified_bilinear
from common.get_Sc_bt import DFfragility
from common.export_outputs import print_outputs
from common.damage_to_loss import damage_to_loss
from spo2ida_based.spo2ida_allTfunction import spo2ida_allT
from spo2ida_based.get_spo2ida_parameters import get_spo2ida_parameters
from common.export_outputs import plot_fragility
from common.conversions import get_spectral_ratios
from common.conversions import from_median_to_mean
pi = 3.141592653589793
plt.close("all")
cd = os.getcwd()

# Define OPTIONS
# -------------------------------------------------------------------------------------------------------
# an_type (Analysis type):
# 0: Elastoplastic Pushover, inelastic displacement obtained from Ruiz-Garcia Miranda [1]
# 1: Any pushover shape (bilinear, trilinear, quadrilinear) inelastic displacement obtained from spo2ida [2]
#
# in_type (input_type):
# 0: input is idealised pushover curve
# 1: input is pushover curve
#
# vulnerability: 1
# 0: stop at fragility curves
# 1: derive vulnerability curves from fragility curves
#
# plotflag (4 integers for each plot: idealised pushover curve, 16%-50%-84% ida curves, fragility curves, vulenarbility curve)
# 1: plot
# 0: don't plot
#
# linew: line width for plots
# fontsize : fontsize used for labels, graphs etc.
# g = units of g compatible with T and displacement
# units: list of 3 strings defining the displacements, forces and Sa units as ['[kN]' '[m]' ['m/s^2]]
# iml = array of intensity measure level used to discretise the fragility curves and get discrete vulnerability curves
#
# Variables for spo2ida procedure:
# pw: pinching model weight (between 0 and 1)
# filletstyle
# 0,1,2 = don't use spline curve to connect different branches of ida curve
# 3: use spline (reccomended)
#
# N: number of points per segment od IDa curve derived with spo2ida
# MC: number of Monte Carlo simulations to account for uncertainty in damage thresholds

an_type = 0
in_type = 1
vulnerability = 1
g = 9.81
iml = np.linspace(0.1,2,50)

plotflag = [1, 0, 1, 1]
linew = 2
fontsize = 10
units = ['[m]', '[kN]', '[g]']
pw = 1
filletstyle = 3
N = 10
MC = 10
Tc = 0.5
Td = 1.8

#-------------------------------------------------------------------------------------------------------
# READ DATA from csv input file and PROCESS DATA:
# First period, first modal participation factor, weights, roof disp. at each limit state, idealised pushover parameters, Dispersions
[T, Gamma, w, dcroof, SPO, bUthd, noBlg] = read_data(in_type,an_type,linew,fontsize,units,plotflag[0])

# Obtain ratios between spectral acceleration at average period and spectral accelerations 
# at the period of each building, in order to produce fragility and vulnerability with a common IM
# and to be able to combine them
Tav = np.average(np.array(T),weights = w) # average period
Sa_ratios = get_spectral_ratios(Tav,T) # ratios

# <codecell>

#-------------------------------------------------------------------------------------------------------
# DERIVE FRAGILITY
plt.close("all")
if vulnerability == 0:
    allSa, allbTSa, allLR50, allbLR = [],[],[],[]     
    for blg in range(0,noBlg):
        # Derive median Sa value (median of Sa) of capacity for each Limit State and corresponding overall dispersion std(log(Sa))
        if an_type==0: # Ruiz-Garcia Miranda's method
            [SaT50, bTSa] = simplified_bilinear(T[blg], Gamma[blg], dcroof[blg], SPO[blg], bUthd[blg], g)
        elif an_type==2: # Dolsek and Fajfar's method
            [mc,a,ac,r,mf] = get_spo2ida_parameters(SPO[blg], T[blg], Gamma[blg]) # Convert MDoF into SDoF
            [SaT50,bTSa] = DFfragility(T[blg], Gamma[blg], dcroof[blg], SPO[blg], bUthd[blg], mc, r, g, Tc, Td)
        else: # Vamvatsikos and Cornell's method
            [mc,a,ac,r,mf] = get_spo2ida_parameters(SPO[blg], T[blg], Gamma[blg]) # Convert MDoF into SDoF
            [idacm, idacr] = spo2ida_allT(mc,a,ac,r,mf,T[blg],pw,plotflag[1],filletstyle,N,linew,fontsize) # apply SPO2IDA procedure
            [SaT50,bTSa] = spo2ida(idacm, idacr, mf, T[blg], Gamma[blg], g, dcroof[blg], SPO[blg], bUthd[blg], MC)
        
        # Converting the Sa(T1) to Sa(Tav), the common IM
        SaTlogmean_av, bTSa_av = np.log(SaT50)*Sa_ratios[blg], np.array(bTSa)*Sa_ratios[blg]
        allSa.append(SaTlogmean_av)
        allbTSa.append(bTSa_av)

    # Combine the fragility of each building in a single lognormal curve with
    # mean = weighted_mean(means) and std = SRSS(weighted_std(means),weighted_mean(stds))
    log_meanSa, log_stSa = [],[]
    for i in range(0,len(dcroof[0])):
        SaLS = [ele[i] for ele in allSa]
        StdSaLS = [ele[i] for ele in allbTSa]
        log_meanSa.append(np.average(SaLS,weights = w)) # weighted log-mean
        log_stSa.append(np.sqrt(np.sum(w*(np.power((SaLS-log_meanSa[i]),2)+np.power(StdSaLS,2))))) # weighted log-std (dispersion)

#-------------------------------------------------------------------------------------------------------
# PLOT AND EXPORT fragility
    # Plot fragility curves
    plt.close("all")
    if plotflag[2]:
        plot_fragility(iml,np.exp(log_meanSa),log_stSa,linew,fontsize,units)   
    # Export fragility parameters (mu and cov of Sa) to csv
    # from log-mean to mean and from dispersion to cofficient of variation
    [meanSa, stSa] = from_median_to_mean(np.exp(log_meanSa),log_stSa)
    cov = np.divide(stSa,meanSa)
    output_path = cd+'/outputs/fragility_parameters.csv'
    header = ['DS', 'mean', 'coefficient of variation']
    n_lines = len(meanSa)
    DS = range(len(meanSa)+1)
    col_data = [DS[1:], meanSa, cov]
    print_outputs(output_path,header,n_lines,col_data)

#-------------------------------------------------------------------------------------------------------
# DERIVE VULNERABILITY curves from damage-to-loss ratios
else: # vulnerability == 1
    LR50s = []
    bLRs = []
    for blg in range(0,noBlg):
        LRs = []
        if an_type != 1:
            if an_type == 0:
                # Ruiz-Garcia Miranda's method the uncertainty in the damage criteria is already included in the total dispersion in a simplified way
                [SaT50, bTSa] = simplified_bilinear(T[blg], Gamma[blg], dcroof[blg], SPO[blg], bUthd[blg], g)
            else:
                [mc,a,ac,r,mf] = get_spo2ida_parameters(SPO[blg], T[blg], Gamma[blg]) # Convert MDoF into SDoF
                [SaT50,bTSa] = DFfragility(T[blg], Gamma[blg], dcroof[blg], SPO[blg], bUthd[blg], mc, r, g, Tc, Td)
            # Conversions
            SaTlogmean_av, bTSa_av = np.log(SaT50)*Sa_ratios[blg], np.array(bTSa)*Sa_ratios[blg]
            LR50 = damage_to_loss(SaTlogmean_av,bTSa_av,iml,linew,fontsize,units)
            LRs.append(LR50)
        else:
            # Vamvatsikos and Cornell's method
            if np.array(bUthd[blg]).any() > 0: 
            # If only some damage criteria have uncertainty
                dc_sample = np.zeros_like(dcroof[blg])
                # Monte Carlo realisations of damage criteria e derivation of fragility curves 
                # for each sample, mean loss ratios are computed at each IML from all samples 
                for i in range (0,MC):
                    for j in range(0,len(dcroof[blg])):
                        if bUthd[blg][j]>0:
                            dc_sample[j] = np.random.lognormal(np.log(dcroof[blg][j]),bUthd[blg][j],1)
                        else:
                            dc_sample[j] = dcroof[blg][j]               
                    [mc,a,ac,r,mf] = get_spo2ida_parameters(SPO[blg], T[blg], Gamma[blg])
                    [idacm, idacr] = spo2ida_allT(mc,a,ac,r,mf,T[blg],pw,[0],filletstyle,N,linew,fontsize)
                    [SaT50, bTSa] = spo2ida(idacm, idacr, mf, T[blg], Gamma[blg], g, dc_sample, SPO[blg], np.zeros_like(bUthd[blg]), MC)
                    # Conversions
                    SaTlogmean_av, bTSa_av = np.log(SaT50)*Sa_ratios[blg], np.array(bTSa)*Sa_ratios[blg]                  
                    LRs.append(damage_to_loss(SaTlogmean_av,bTSa_av,iml,linew,fontsize,units))
            else: 
            # If any damage criteria have uncertainty
                [mc,a,ac,r,mf] = get_spo2ida_parameters(SPO[blg], T[blg], Gamma[blg])
                [idacm, idacr] = spo2ida_allT(mc,a,ac,r,mf,T[blg],pw,[0],filletstyle,N,linew,fontsize)
                [SaT50, bTSa] = spo2ida(idacm, idacr, mf, T[blg], Gamma[blg], g, dcroof[blg], SPO[blg], bUthd[blg], MC)
                # Conversions
                SaTlogmean_av, bTSa_av = np.log(SaT50)*Sa_ratios[blg], np.array(bTSa)*Sa_ratios[blg] 
                # Reconvert to median and dispersion
                LRs.append(damage_to_loss(SaTlogmean_av,bTSa_av,iml,linew,fontsize,units))
        # Define Vulnerability curve for each building
        LR50s.append(np.mean(LRs,axis = 0))
        bLRs.append(np.std(LRs, axis = 0))
        
    # Combine the Loss Ratios of each building in a single LR for each iml, 
    # with mean = mean(LRs) and std = std(LRs)
    lr = np.array(LR50s)
    bLRs = np.array(bLRs)
    LR50 = np.average(lr, axis=0, weights=w) # weighted mean
    w2 = np.array([np.repeat(ele,len(iml)) for ele in w])
    bLR = np.sqrt(np.sum(np.multiply(np.power(bLRs,2)+np.power(lr-LR50,2),w2),axis=0))

#-------------------------------------------------------------------------------------------------------
# PLOT AND EXPORT vulnerability   
    # Plot mean Vulnerability curve
    if plotflag[3]:
        cd = os.getcwd()
        plt.close("all")
        plt.plot(iml,LR50, marker='o', linestyle='None',label = 'LR')
        plt.xlabel('Spectral acceleration at T elastic, Sa(Tel) '+units[2],fontsize = fontsize)
        plt.ylabel('Loss Ratio',fontsize = fontsize)
        plt.suptitle('Vulnerability Curve',fontsize = fontsize)
        plt.legend(loc='lower right',frameon = False)
        plt.savefig(cd+'/outputs/vulnerability_curve.png')
        plt.show()

    # Export discrete vulnerability curve to csv
    output_path = cd+'/outputs/vulnerability_curve.csv'
    header = ['Intensity Measure Level','mean Loss Ratio','coefficient of variation Loss Ratio']
    n_lines = len(LR50)
    col_data = [iml, LR50, bLR/LR50] # mean and coefficient of variation of LR
    print_outputs(output_path,header,n_lines,col_data)



