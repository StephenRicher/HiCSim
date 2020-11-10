#!/usr/bin/env python3

""" Write distribution and frequency of TUs within polymer """

import sys
import argparse
import pandas as pd
from typing import List
from utilities import setDefaults, readJSON

__version__ = '1.0.0'


def writeTUdistribution(infiles: List, out: str):

    allBeadDistribution = {}
    for i, file in enumerate(infiles):
        beads = readJSON(file)
        # Set all beads in first iteration.
        if i == 0:
            for atom in beads['DNA']:
                allBeadDistribution[atom] = 0
        for TU in beads['TU']:
            allBeadDistribution[TU] += 1

    allBeadDistribution = pd.DataFrame(allBeadDistribution, index=[0])
    allBeadDistribution.to_csv(out, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infiles', nargs='+', help='Atom JSON file.')
    parser.add_argument('--out', default=sys.stdout, help='Output filename')
    
    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(writeTUdistribution(**vars(args)))
