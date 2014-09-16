# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 16:03:26 2014

@author: chiaracasotto
"""
import csv

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