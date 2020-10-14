#!/usr/bin/env python3

""" Determine TU activation state at each time point """

import sys
import gzip
import json
import logging
import argparse
import fileinput
import numpy as np
import pandas as pd
from utilities import readCustom
from scipy.spatial.distance import cdist


__version__ = '1.0.0'

def main(file: str, atomGroups: str, distance: float, timestep: float, **kwargs) -> None:

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)

    with gzip.open(file, 'rt') as fh1, gzip.open(file, 'rt') as fh2:
        allDistances = []
        while True:
            try:
                TUs = readCustom(fh1, includeIDs=atomGroupsDict['TU'])
                # Read TFs and exclude type 3 (inactive TF)
                TFas = readCustom(fh2, includeIDs=atomGroupsDict['TF'], excludeTypes=['3'])
                # Compute distance of active TFs to each TU bead
                result = cdist(TUs['atoms'], TFas['atoms'], 'euclidean')
                # Find closest monomer to each bead
                minDistance = np.amin(result, axis=1) < distance
                allDistances.append(minDistance)
            except EOFError:
                break
        # Convert list of distance arrays to 2D numpy matrix
        allDistances = np.vstack(tuple(allDistances))

    TU_indexes = atomGroupsDict['TU']

    # Convert to pandas so we can label columns with TU atom indexes
    allDistances = pd.DataFrame(data=allDistances, columns=TU_indexes)

    # Convert to long format
    allDistances = allDistances.stack().reset_index()
    allDistances.columns = ['timepoint', 'TU', 'activated']
    allDistances['timepoint'] *= timestep
    allDistances.to_csv(sys.stdout, index=False)


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
        'file',
        help='Input gzipped custom lammps simulation file.')
    custom.add_argument(
        '--distance', default=1.8, type=float,
        help='Distance threshold for active transcription (default: %(default)s)')
    custom.add_argument(
        '--timestep', default=1, type=int,
        help='Number of timesteps per sampling (default: %(default)s)')
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
