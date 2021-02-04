#!/usr/bin/env python3

""" Read all BEAD files and map type to an ID number """

import sys
import json
import logging
import argparse
import fileinput
from typing import List
from argUtils import setDefaults, createMainParent


__version__ = '1.0.0'


def extractAtomTypes(beads: List):
    beadIDs = {'TFa': 1, 'TFi': 2, 'N': 3}
    id = 4 # Start at 4 since 1, 2, 3 are reserved defined types (TFa, TFi, N)
    with fileinput.input(beads) as fh:
        for line in fh:
            line = line.strip()
            if line in ['TFa', 'TFi']:
                logging.error(
                    'TFa and Tfi beads now allowed in polymer bead files.')
                return 1
            if line not in beadIDs:
                beadIDs[line] = id
                id += 1
    json.dump(beadIDs, sys.stdout)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=extractAtomTypes)
    parser.add_argument(
        'beads', nargs='*', help='Input bead files (defautlt: stdin)')
    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
