#!/usr/bin/env python3

import re
import sys
import json
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

def setDefaults(parser):
    """ Set logging config and return args with associated function """

    args = parser.parse_args()
    logFormat='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    try:
        logging.basicConfig(level=args.verbose, format=logFormat)
        del args.verbose
    except AttributeError:
        logging.basicConfig(format=logFormat)
        pass
    try:
        function = args.function
        del args.function
    except AttributeError:
        parser.print_help()
        sys.exit()

    return args, function


def createMainParent(verbose=True, version=None):
    """ Create parser of verbose/version to be added to parser/subparsers. """
    parent = argparse.ArgumentParser(add_help=False)
    if version:
        parent.add_argument('--version', action='version',
            version=f'%(prog)s {version}')
    if verbose:
        parent.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')
    return parent


def transformScore(score, transform='none'):
    """ Apply specified transformation """

    assert transform in ['none', 'sqrt', 'log']
    if transform == 'sqrt':
        return np.sqrt(score)
    elif transform == 'log':
        return np.log(score)
    else:
        return score


def bedHeader(line):
    """ Return True for empty lines of header strings """

    line = line.strip()
    if not line or line.startswith(('browser', 'track', '#')):
        return True
    else:
        return False

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


def getBead(pos, nbases):
    """ Return bead corresponding to nbases """

    return pos // nbases

def pearsonr_pval(x,y):
    try:
        return pearsonr(x,y)[1]
    except ValueError:
        return np.nan


def countPair(x, y):
    """ Return count of valid pairs (both not nan) """

    # Indices where both x and y are NOT np.nan
    validIndices = np.intersect1d(
        np.where(~np.isnan(x)),
        np.where(~np.isnan(y)))
    return len(validIndices)
