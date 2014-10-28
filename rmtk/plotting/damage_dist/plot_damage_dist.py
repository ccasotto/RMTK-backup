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

def plot_damage_dist(damage_file, taxonomy_list=[], log_scale_x=True, log_scale_y=True, export_png=False):
    '''
    Plots the damage distribution for the specified taxonomies
    '''
    metadata, asset_refs, damage_dist = parsedd.parse_damage_file(damage_file)
    if taxonomy_list:
        asset_refs = taxonomy_list
    for ref in asset_refs:
        loss, poe = damage_dist[ref]
        fig = pyplot.figure(figsize = (16, 9))
        if log_scale_x:
            if log_scale_y:
                pyplot.loglog(loss, poe, '-r', linewidth = 2, label = 'Asset ' + ref)
            else:
                pyplot.semilogx(loss, poe, '-r', linewidth = 2, label = 'Asset ' + ref)
        elif log_scale_y:
            pyplot.semilogy(loss, poe, '-r', linewidth = 2, label = 'Asset ' + ref)
        else:
            pyplot.plot(loss, poe, '-r', linewidth = 2, label = ref)
        pyplot.title('Loss curve (' + metadata['lossType'] + ' losses)', fontsize = 20)
        pyplot.legend(loc = "upper right", bbox_to_anchor = (1,1))
        pyplot.xlabel('Loss (' + metadata['unit'] + ')', fontsize = 16)
        pyplot.ylabel('Probability of exceedance in ' + metadata['investigationTime'] + ' years', fontsize = 16)
        pyplot.grid(b=None, which='major', axis='both', color='LightGray', linestyle='--', linewidth=.5)
        if export_png:
            pyplot.savefig(ref, format = 'png')

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