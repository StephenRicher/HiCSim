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
from itertools import combinations
from scipy.spatial.distance import cdist, pdist


__version__ = '1.0.0'

def main(file: str, atomGroups: str, outDistances: str, distance: float, timestep: float, **kwargs) -> None:

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)
    names1, names2 = generatePairwiseNames(atomGroupsDict['TU'])

    with gzip.open(file, 'rt') as fh1, gzip.open(file, 'rt') as fh2:
        allActiveTUs = []
        allTUpairs = []
        while True:
            try:
                TUs = readCustom(fh1, includeIDs=atomGroupsDict['TU'])
                # Read TFs and exclude type 3 (inactive TF)
                TFas = readCustom(fh2, includeIDs=atomGroupsDict['TF'], excludeTypes=['3'])
                # Compute distance of active TFs to each TU bead
                result = cdist(TUs['atoms'], TFas['atoms'], 'euclidean')
                # Find closest monomer to each bead
                activeTUs = np.amin(result, axis=1) < distance
                allActiveTUs.append(activeTUs)

                activeTUpos = []
                inactiveTUpos = []
                for i, TU in enumerate(TUs['atoms']):
                    if activeTUs[i]:
                        activeTUpos.append(TU)
                        inactiveTUpos.append([np.nan, np.nan, np.nan])
                    else:
                        activeTUpos.append([np.nan, np.nan, np.nan])
                        inactiveTUpos.append(TU)

                TUpairs = pd.DataFrame(
                    {'TU1'      : names1                         ,
                     'TU2'      : names2                         ,
                     'active'   : pdist(activeTUpos,   'euclidean'),
                     'inactive' : pdist(inactiveTUpos, 'euclidean')})
                allTUpairs.append(TUpairs)
            except EOFError:
                break

    # Average active TU-TU pair distances across timepoints and save
    allTUpairs = pd.concat(allTUpairs)
    allTUpairs = allTUpairs.groupby(['TU1', 'TU2']).mean().reset_index()
    allTUpairs.to_csv(outDistances, index=False)

    # Convert to pandas so we can label columns with TU atom indexes
    allActiveTUs = pd.DataFrame(
        allActiveTUs, columns=atomGroupsDict['TU']).stack().reset_index()
    allActiveTUs.columns = ['timepoint', 'TU', 'activated']
    allActiveTUs['timepoint'] *= timestep
    allActiveTUs.to_csv(sys.stdout, index=False)


def readJSON(file):
    """ Read JSON encoded data to dictionary """
    with open(file) as fh:
        return json.load(fh)


def generatePairwiseNames(names):
    names1 = []
    names2 = []
    names = combinations(names, 2)
    for name1, name2 in names:
        names1.append(name1)
        names2.append(name2)
    return names1, names2


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
        '--outDistances', default='TU-activePairDistance.csv',
        help='File to write active TU pair distances.')
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
