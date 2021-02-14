#!/usr/bin/env python3

""" Script to read LAMMPS simulation and generate contact frequency matrix """


import sys
import logging
import argparse
import numpy as np
import pandas as pd
from typing import List
from scipy.sparse import save_npz, csc_matrix
from scipy.spatial.distance import pdist, squareform
from argUtils import positiveInt
from utilities import getAtomCount, readJSON, pdistPeriodic, setDefaults


__version__ = '1.0.0'


def createContactMatrix(
        file: str, atomGroups: str, periodic: bool, seed: int,
        dimensions : List, out: str, distance: float):
    np.random.seed(seed)
    if periodic and not dimensions:
        logging.warning(
            '"--periodic" set but ignored as "--dimensions" not set.')
    elif dimensions and not periodic:
        logging.warning(
            '"--dimensions" set but ignored as "--periodic" not set.')

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)
    atomCount = getAtomCount(atomGroupsDict)

    # Read simulation file 1 timestep at a time and filter as needed
    fullSim = pd.read_csv(file, chunksize=atomCount)

    sqdistance = distance**2
    contacts = 0
    for timestep in fullSim:
        timestep = timestep[timestep['ID'].isin(atomGroupsDict['DNA'])]
        xyz = timestep[['x', 'y', 'z']].to_numpy()
        if periodic and dimensions:
            distance = pdistPeriodic(xyz, dimensions, sqeuclidean=True)
        else:
            distance = pdist(xyz, 'sqeuclidean')

        # Transform to probabilities
        distance = np.exp(-distance / sqdistance)
        # Accept contact probabilities
        distance = distance > np.random.rand(*distance.shape)
        contacts += distance


    save_npz(out, csc_matrix(squareform(contacts)))


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'atomGroups',
        help='Atom group assignments in JSON format.')
    parser.add_argument(
        'file', nargs='?', default=[],
        help='Input custom lammps simulation file (default: stdin)')
    parser.add_argument(
        '--periodic', action='store_true',
        help='Compute euclidean distance considering 3D periodic boundary '
             'condition. Ignored if "--dimensions" not set.')
    parser.add_argument(
        '--dimensions', metavar=('X', 'Y', 'Z'), nargs=3, type=float,
        help='Simulation dimensions. Ignored if "--periodic" not set.')
    parser.add_argument(
        '--distance', default=3, type=float,
        help='Distance threshold for determing contact probability '
             '(default: %(default)s)')
    parser.add_argument(
        '--seed', default=None, type=positiveInt,
        help='Non-negative integer for seeding random placement '
             'positions (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True, help='Output contact NPZ matrix.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(createContactMatrix(**vars(args)))
