# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 16:03:26 2014

@author: chiaracasotto
"""
import csv
import os
import matplotlib.pyplot as plt
import scipy.stats as stat
import numpy as np

def print_outputs(output_file,header,n_lines,col_data):

    # Open and read to write the file a just created
    outfile = open(output_file, 'wt')
    writer = csv.writer(outfile, delimiter=',')

    # Writes the header
    writer.writerow(header)
    # Loop over the data vector
    for j in range(0, n_lines):
        dat = [ele[j] for ele in col_data]
        #dat = [DS, np.log(Sa50[j]), bTSa[j]] 
        # Writes each line 
        writer.writerow(dat)
        # Close the file 
    outfile.close()

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