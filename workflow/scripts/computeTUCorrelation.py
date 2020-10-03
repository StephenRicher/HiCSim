#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency matrix """

import sys
import json
import logging
import argparse
import fileinput
import itertools
import numpy as np
import pandas as pd
from utilities import readCustom
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr


__version__ = '1.0.0'

def main(file: str, atomGroups: str, out: str, distance: float, **kwargs) -> None:

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)

    with fileinput.input(file) as fh:
        allDistances = []
        while True:
            try:
                TUs = readCustom(fh, includeIDs=atomGroupsDict['TU'])
                # Read TFs and exclude type 3 (inactive TF)
                TFas = readCustom(fh, includeIDs=atomGroupsDict['TF'], excludeTypes=['3'])
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

    # Compute Pearson correlation coefficient and convert to long format
    correlations = allDistances.corr(method='pearson')
    correlations = correlations.stack().reset_index()
    correlations.columns = ['row', 'column', 'r']

    # Combpute corresponding p-values and convert to seperate dataframe
    pvalues = allDistances.corr(method=pearsonr_pval)
    pvalues = pvalues.stack().reset_index()
    pvalues.columns = ['row', 'column', 'p-value']

    # Combine to include p-values
    allDistances = pd.merge(correlations, pvalues)

    allDistances.to_csv(out, index=False)


def pearsonr_pval(x,y):
    return pearsonr(x,y)[1]


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
        '--distance', default=1.8, type=float,
        help='Distance threshold for active transcription (default: %(default)s)')
    custom.add_argument(
        '--out', default='TUcorrelation.csv',
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
