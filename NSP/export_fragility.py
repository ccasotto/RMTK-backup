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

def export_fragility(vuln, plot_feature, x, y,leg):
    plotflag, linew, fontsize, units, iml = plot_feature[0:5]
    cd = os.getcwd()
    if vuln == 0:
        if plotflag[2]:
            plot_fragility(iml,np.exp(x),y,linew,fontsize,units,leg)
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

def plot_fragility(iml,Sa50,bTSa,linew,fontsize,units,leg):
    # INPUT: Sa50 is the median of iml, while bTSa is the dispersion, that is to say the std(log(iml))
    colours = ['b','r','g','k','c','y']
    cd = os.getcwd()
    txt = []
    for q in range(0,len(Sa50)):
        txt.append('Damage State '+str(q+1))
        y = stat.norm(np.log(Sa50[q]),bTSa[q]).cdf(np.log(iml))
        if leg == 'off':
            plt.plot(iml,y,'--',color=colours[q],linewidth=linew)
        else:
            plt.plot(iml,y,color=colours[q],linewidth=linew)

    plt.hold(True)
    if leg == 'on':
        plt.legend(txt,loc='lower right',frameon = False,fontsize = fontsize)
        plt.xlabel('Spectral acceleration at T elastic, Sa(Tel) [g]',fontsize = fontsize)
        plt.ylabel('Probabilty of Exceedance',fontsize = fontsize)
        plt.suptitle('Fragility Curves',fontsize = fontsize)
        plt.savefig(cd+'/outputs/fragility_curves.png')
        plt.show()
