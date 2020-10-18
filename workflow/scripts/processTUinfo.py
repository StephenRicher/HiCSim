#!/usr/bin/env python3

""" Determine TU activation state at each time point """

import sys
import json
import logging
import argparse
import fileinput
import numpy as np
import pandas as pd
from utilities import readCustom2
from itertools import combinations
from scipy.spatial.distance import cdist, pdist


__version__ = '1.0.0'

def main(file: str, atomGroups: str, distance: float, outTU: str, outPolymer: str, outPairDistance: str, **kwargs) -> None:

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)
    names1, names2 = generatePairwiseNames(atomGroupsDict['TU'])

    with fileinput.input(file) as fh:
        allTUinfo = []
        allTUpairs = []
        i = 0
        while True:
            try:
                # Read simulation data for a timepoint
                sim = readCustom2(fh)

                # Extract data of only TUs
                TUdata = sim.loc[sim['id'].isin(atomGroupsDict['TU'])]
                TUs = TUdata[['xPos', 'yPos', 'zPos']].to_numpy()

                # Read TFs and exclude type 3 (inactive TF)
                TFas = sim.loc[(sim['id'].isin(atomGroupsDict['TF']))
                    & (sim['type'] != '3'),
                    ['xPos', 'yPos', 'zPos']].to_numpy()

                # Compute distance of active TFs to each TU bead
                result = cdist(TUs, TFas, 'euclidean')

                TUdata = TUdata.copy()
                TUdata['TFnear'] = np.amin(result, axis=1)
                TUdata['active'] = TUdata['TFnear'] < distance
                allTUinfo.append(TUdata)

                # Write all DNA position, for each timepoint, to file
                # Appended to file to save memory
                if outPolymer:
                    header = True if i == 0 else False
                    PolymeInfo = sim.loc[sim['id'].isin(atomGroupsDict['DNA'])]
                    PolymeInfo.to_csv(
                        outPolymer, mode='a', header=header, index=False)

                activeTUpos = []
                inactiveTUpos = []
                activeTUs = TUdata['active'].to_list()
                for i, TU in enumerate(TUs):
                    if activeTUs[i]:
                        activeTUpos.append(TU)
                        inactiveTUpos.append([np.nan, np.nan, np.nan])
                    else:
                        activeTUpos.append([np.nan, np.nan, np.nan])
                        inactiveTUpos.append(TU)

                allTUpairs.append(pd.DataFrame(
                    {'TU1'      : names1                         ,
                     'TU2'      : names2                         ,
                     'active'   : pdist(activeTUpos,   'euclidean'),
                     'inactive' : pdist(inactiveTUpos, 'euclidean')}))

                i += 1
            except EOFError:
                break

    # Average active TU-TU pair distances across timepoints and save
    allTUpairs = pd.concat(allTUpairs)
    allTUpairs = allTUpairs.groupby(['TU1', 'TU2']).mean().reset_index()
    allTUpairs.to_csv(outPairDistance, index=False)

    pd.concat(allTUinfo).to_csv(outTU, index=False)


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
        'file', nargs='?', default=[],
        help='Input custom lammps simulation file (default: stdin)')
    custom.add_argument(
        '--distance', default=1.8, type=float,
        help='Distance threshold for active transcription (default: %(default)s)')
    custom.add_argument(
        '--outPairDistance', default='TU-pairDistance.csv.gz',
        help='File to write active TU pair distances (default: %(default)s)')
    custom.add_argument(
        '--outTU', default='TU-info.csv.gz', help='File to TU info.')
    custom.add_argument(
        '--outPolymer', default='Polymer-info.csv.gz',
        help='File to polymer info.')
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
