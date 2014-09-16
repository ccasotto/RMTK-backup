# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 16:43:42 2014

@author: chiaracasotto
"""
import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import os
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