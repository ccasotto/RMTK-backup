'''
Post-process damage calculation data to extract the damage distribution
for different taxonomies
'''

import os
import csv
import argparse
import numpy as np
from lxml import etree
from collections import OrderedDict

xmlNRML='{http://openquake.org/xmlns/nrml/0.4}'
xmlGML = '{http://www.opengis.net/gml}'

def parse_single_damage_dist(element):
    '''
    Reads the dmgDist element to return the longitude, latitude and 
    poes and losses
    '''
    for e in element.iter():
        ref = element.attrib.get('assetRef')
        if e.tag == '%spos' % xmlGML:
            coords = str(e.text).split()
            lon = float(coords[0])
            lat = float(coords[1])
        elif e.tag == '%spoEs' % xmlNRML:
            poes = str(e.text).split()
            poes = map(float, poes)
        elif e.tag == '%slosses' % xmlNRML:
            losses = str(e.text).split()
            losses = map(float, losses)
        else:
            continue
    return lon, lat, ref, poes, losses


def parse_damage_states(element):
    '''
    Returns the statistics
    '''
    damage_states = {}
    damage_states = element.attrib.get('damageStates')
    return damage_states


def parse_damage_file(input_file):
    '''
    Reads an xml damage dist file and returns a dictionary with
    taxonomies as keys and (mean[ds], stddev[ds]) as values
    '''
    damage_states = {}
    taxonomies = []
    for _, element in etree.iterparse(input_file):
        if element.tag == '%sdamageStates' % xmlNRML:
            damage_states = str(element.text).split()
        elif element.tag == '%slossCurve' % xmlNRML:
            lon, lat, ref, poe, loss = parse_single_damage_dist(element)
            asset_refs.append(ref)
            damage_dists[ref] = loss, poe
        else:
            continue
    return damage_states, asset_refs, damage_dists


def set_up_arg_parser():
    """
    Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert NRML damage distribution files'
            'to csv files. To run just type: python parse_damage_dists.py '
            '--input-file=PATH_TO_DAMAGE_DIST_NRML_FILE ', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--input-file',
        help='path to damage distribution NRML file (Required)',
        default=None,
        required=True)

    return parser

if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        damage_states, taxonomies, damage_dists = parse_damage_file(args.input_file)