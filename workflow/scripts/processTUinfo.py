#!/usr/bin/env python3

""" Determine TU activation state at each time point """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
from utilities import getAtomCount, readJSON
from scipy.spatial.distance import cdist


__version__ = '1.0.0'

def main(file: str, atomGroups: str, distance: float, out: str, **kwargs) -> None:

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)
    atomCount = getAtomCount(atomGroupsDict)

    # Read simulation file 1 timestep at a time
    info = pd.read_csv(file,
        usecols=['time', 'id', 'type', 'x', 'y', 'z'], chunksize=atomCount)

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
        '--out', default=sys.stdout, help='File to TU info (default: stdout)')
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
