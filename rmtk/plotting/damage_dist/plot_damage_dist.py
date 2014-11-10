'''
Post-process damage calculation outputs to plot damage distibution charts
'''

import os
import csv
import argparse
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot
import parse_damage_dist as parsedd

xmlNRML = '{http://openquake.org/xmlns/nrml/0.4}'
xmlGML = '{http://www.opengis.net/gml}'

def parse_taxonomy_file(taxonomy_file):
    '''
    Reads a txt file with a list of taxonomies and returns an array
    '''
    taxonomy_list = []
    file = open(assets_file,'r')
    taxonomy_list = file.read().split(',')
    file.close()
    return taxonomy_list

def plot_damage_dist(damage_file, taxonomy_list=[], plot_3d=False, export_png=False):
    '''
    Plots the damage distribution for the specified taxonomies
    '''
    taxonomies, damage_states, damage_dist_tax = parsedd.parse_damage_file(damage_file)
    print damage_states
    if taxonomy_list:
        taxonomies = taxonomy_list
    for tax in taxonomies:
        means = []
        stddevs = []
        damage_dist = damage_dist_tax[tax]
        for ds in damage_states:
            dd = damage_dist[ds]
            means.append(dd[0])
            stddevs.append(dd[1])
        fig = pyplot.figure(figsize = (16, 9))

        N = len(damage_states)
        ind = np.arange(N)  # the x locations for the groups
        error_config = {'ecolor': '0.3', 'linewidth': '2'}
        bar_width = 0.3
        padding_left = 0
        pyplot.bar(ind+padding_left, height=means, width=bar_width, yerr=stddevs, error_kw=error_config, color='IndianRed', linewidth=1.5)
        pyplot.title('Damage distribution (' + tax + ')', fontsize = 20)
        pyplot.xlabel('Damage state', fontsize = 16)
        pyplot.ylabel('Number of assets in damage state', fontsize = 16)
        pyplot.xticks(ind+padding_left+bar_width/2., damage_states)
        pyplot.margins(.25,0)
        if export_png:
            pyplot.savefig(tax, format = 'png')
        pyplot.show()
        pyplot.clf()

def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description = 'Convert NRML loss curves file to tab delimited '
            ' .txt files. Inside the specified output directory, create a .txt '
            'file for each stochastic event set.'
            'To run type: python plot_damage_dist.py '
            '--input-file = PATH_TO_LOSS_CURVE_NRML_FILE '
            '--taxonomy_list = LIST_OF_ASSETS ', add_help = False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action = 'help')
    flags.add_argument('--input-file',
        help = 'path to loss curves NRML file (Required)',
        default = None,
        required = True)
    flags.add_argument('--assets-list',
        help = 'path to loss curves NRML file (Required)',
        default = [],
        required = False)
    return parser

if __name__ ==  "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        if args.taxonomy_list:
            plot_damage_dist(args.input_file, taxonomy_list, log_scale_x=True, log_scale_y=True, export_png=False)
        else:
            plot_damage_dist(args.input_file, taxonomy_list, log_scale_x=True, log_scale_y=True, export_png=False)
    else:
        parser.print_usage()