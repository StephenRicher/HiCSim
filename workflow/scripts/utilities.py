#!/usr/bin/env python3

import sys
import re
import json
import logging
import argparse
import pandas as pd
import numpy as np


def setDefaults(parser, verbose=True, version=None):
    """ Add version and verbose arguments to parser """

    if version:
        parser.add_argument('--version', action='version',
            version=f'%(prog)s {version}')
    if verbose:
        parser.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')

        args = parser.parse_args()
        logFormat='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
        logging.basicConfig(level=args.verbose, format=logFormat)
        del args.verbose
    else:
        args = parser.parse_args()

    return args


#https://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co/11109336#11109336
def cdistPeriodic(x0, x1, dimensions, sqeuclidean=False):
    delta = np.abs(x0[:, np.newaxis] - x1)
    dimensions = np.array(dimensions)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    sqdistances = (delta ** 2).sum(axis=-1)
    if sqeuclidean:
        return sqdistances
    else:
        return np.sqrt(sqdistances)


def pdistPeriodic(x, dimensions, sqeuclidean=False):
    # Retrive pairwise indices
    r,c = np.triu_indices(len(x),1)
    # Subtract only non-repeating pairwise
    delta = np.abs(x[r] - x[c])
    dimensions = np.array(dimensions)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    sqdistances = (delta ** 2).sum(axis=-1)
    if sqeuclidean:
        return sqdistances
    else:
        return np.sqrt(sqdistances)


def getAtomCount(atomGroups, key=None):
    """ Read dictionary of atom groups and return total number of atoms """

    uniqueAtoms = set()
    if key:
        atomGroups = {key : atomGroups[key]}
    for atoms in atomGroups.values():
        uniqueAtoms.update(atoms)
    return len(uniqueAtoms)


def readJSON(file):
    """ Read JSON encoded data to dictionary """
    with open(file) as fh:
        return json.load(fh)


def npz(value):
    ''' Ensure output file ends with '.npz'. '''

    if not value.endswith('.npz'):
        raise argparse.ArgumentTypeError(
            f'Out file {value} must end with: \'.npz\'.')
    else:
        return value


def scale(value):
    ''' Ensure input is numeric or 'None'. '''

    if value == 'None':
        return None
    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f'must be float value or None')


def region(value):
    ''' Validate input for region argument of create_contact_matrix '''

    pattern = re.compile('^[0-9]+-[0-9]+$')
    if not pattern.match(value):
        raise argparse.ArgumentTypeError(
            f'Expected format is START-END e.g 1-100000000')
    regions = {}
    regions['start'] = int(value.split('-')[0])
    regions['end'] = int(value.split('-')[1])

    if not regions['start'] < regions['end']:
        raise argparse.ArgumentTypeError(
            f'Start coordinate {regions["start"]} not less '
            f'than end coordinate {regions["end"]}.')
    else:
        return regions

def coordinates(value):
    ''' Validate input for genomic coordinates  '''

    pattern = re.compile('^[^:-]+:[0-9]+-[0-9]+$')
    if not pattern.match(value):
        raise argparse.ArgumentTypeError(
            'Expected format is CHR:START-END e.g chr1:1-1000. '
            'Chromosome name cannot contain ": -" .')
    coords = {}
    coords['chr'], coords['start'], coords['end'] = re.sub('[:-]', ' ', value).split()
    coords['start'] = int(coords['start'])
    coords['end'] = int(coords['end'])
    if not coords['start'] < coords['end']:
        raise argparse.ArgumentTypeError(
            f'Start coordinate {coords["start"]} not less '
            f'than end coordinate {coords["end"]}.')
    else:
        return coords


def commaPair(value):
    ''' Split command seperated pair and return as dictionary '''

    value = value.split(',')
    v1 = value[0].strip()
    v2 = value[1].strip()
    if len(v2) != 1:
        raise argparse.ArgumentTypeError(
            f'Masking character {v2} shoud be single character.')
    else:
        return (v1, v2)


def coeff(value):
    ivalue = float(value)
    if not -1 <= ivalue <= 1:
        raise argparse.ArgumentTypeError(f'{value} must be between -1 and 1.')
    return ivalue


def getBead(pos, nbases):
    """ Return bead corresponding to nbases """

    return pos // nbases
