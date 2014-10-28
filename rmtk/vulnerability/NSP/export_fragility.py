# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 16:03:26 2014

@author: chiaracasotto
"""
import os
import matplotlib.pyplot as plt
import scipy.stats as stat
import numpy as np
from common.conversions import from_median_to_mean
from common.print_csv import print_outputs

def export_fragility(vuln, plot_feature, x, y):
    plotflag, linew, fontsize, units, iml = plot_feature[0:5]
    cd = os.getcwd()
    if vuln == 0:
        if plotflag[2]:
            plot_fragility(iml,np.exp(x),y,linew,fontsize,units)   
        # Export fragility parameters (mu and cov of Sa) to csv
        # from log-mean to mean and from dispersion to cofficient of variation
        [meanSa, stSa] = from_median_to_mean(np.exp(x),y)
        cov = np.divide(stSa,meanSa)
        output_path = cd+'/outputs/fragility_parameters.csv'
        header = ['DS', 'mean', 'coefficient of variation']
        n_lines = len(meanSa)
        DS = range(len(meanSa)+1)
        col_data = [DS[1:], meanSa, cov]
        print_outputs(output_path,header,n_lines,col_data)

def plot_fragility(iml,Sa50,bTSa,linew,fontsize,units):
    # INPUT: Sa50 is the median of iml, while bTSa is the dispersion, that is to say the std(log(iml))
    colours = ['b','r','g','k','c','y']
    cd = os.getcwd()
    #texto = ['yielding','collapse','mod']
    for q in range(0,len(Sa50)):
        txt = 'Damage State '+str(q+1)
        y = stat.norm(np.log(Sa50[q]),bTSa[q]).cdf(np.log(iml))
        plt.plot(iml,y,color=colours[q],linewidth=linew,label = txt)
    
    plt.xlabel('Spectral acceleration at T elastic, Sa(Tel) '+units[2],fontsize = fontsize)
    plt.ylabel('Probabilty of Exceedance',fontsize = fontsize)
    plt.suptitle('Fragility Curves',fontsize = fontsize)
    plt.legend(loc='lower right',frameon = False)
    plt.savefig(cd+'/outputs/fragility_curves.png')
    plt.show()