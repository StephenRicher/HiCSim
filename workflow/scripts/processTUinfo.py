#!/usr/bin/env python3

""" Determine TU activation state at each time point """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
from utilities import setDefaults, getAtomCount, readJSON
from scipy.spatial.distance import cdist

__version__ = '1.0.0'


def processTUinfo(file: str, atomGroups: str, distance: float, out: str):

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)
    if ('TU' not in atomGroupsDict) or ('TF' not in atomGroupsDict):
        logging.error('No TU-TF pairs to process in simulation.')
        return 1

    atomCount = getAtomCount(atomGroupsDict)

    # Read simulation file 1 timestep at a time
    info = pd.read_csv(
        file, chunksize=atomCount,
        usecols=['time', 'id', 'type', 'x', 'y', 'z'])

    allTUinfo = []
    for sim in info:
        # Extract data of only TUs
        TUdata = sim.loc[sim['id'].isin(atomGroupsDict['TU'])]
        TUs = TUdata[['x', 'y', 'z']].to_numpy()

        # Read TFs and exclude type 3 (inactive TF)
        TFas = sim.loc[(sim['id'].isin(atomGroupsDict['TF']))
                       & (sim['type'] != '3'), ['x', 'y', 'z']].to_numpy()

        # Compute distance of active TFs to each TU bead
        result = cdist(TUs, TFas, 'euclidean')

        TUdata = TUdata.copy()
        TUdata['TFnear'] = np.amin(result, axis=1)
        TUdata['active'] = TUdata['TFnear'] < distance
        allTUinfo.append(TUdata)

    pd.concat(allTUinfo).to_csv(out, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'atomGroups', help='Atom group assignments in JSON format.')
    parser.add_argument(
        'file', nargs='?', default=[],
        help='Input custom lammps simulation file (default: stdin)')
    parser.add_argument(
        '--distance', default=1.8, type=float,
        help='Distance threshold for transcription (default: %(default)s)')
    parser.add_argument(
        '--out', default=sys.stdout, help='File to TU info (default: stdout)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(processTUinfo(**vars(args)))
