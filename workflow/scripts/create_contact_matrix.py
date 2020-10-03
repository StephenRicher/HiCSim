#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency matrix """

import sys
import json
import logging
import argparse
import fileinput
import numpy as np
from utilities import readCustom
from timeit import default_timer as timer
from scipy.sparse import save_npz, csc_matrix
from scipy.spatial.distance import pdist, squareform

__version__ = '1.0.0'

def main(file: str,  atomGroups: str, outdata: str, distance: float, **kwargs) -> None:

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)

    sqdistance = distance**2
    contacts = 0
    with fileinput.input(file) as fh:
        while True:
            try:
                xyz = readCustom(fh, includeIDs=atomGroupsDict['DNA'])
                contacts += pdist(xyz['atoms'], 'sqeuclidean') < sqdistance
            except EOFError:
                break

    save_npz(outdata, csc_matrix(squareform(contacts)))


def readJSON(file):
    """ Read JSON encoded data to dictionary """
    with open(file) as fh:
        return json.load(fh)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'atomGroups',
        help='Atom group assignments in JSON format.')
    custom.add_argument(
        'file', nargs='?', default=[],
        help='Input custom lammps simulation file (default: stdin)')
    custom.add_argument(
        '--distance', default=3, type=float,
        help='Max contact distance between particles (default: %(default)s)')
    custom.add_argument(
        '--outdata', default='contacts.npz',
        help='Contact matrix output (default: %(default)s)')
    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')
    base.add_argument(
        '--verbose', action='store_const', const=logging.DEBUG,
        default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom])
    args = parser.parse_args()

    log_format='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    logging.basicConfig(level=args.verbose, format=log_format)

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
