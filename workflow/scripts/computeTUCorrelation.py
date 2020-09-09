#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency matrix """

import sys
import json
import logging
import argparse
import itertools
import numpy as np
import pandas as pd
from utilities import read_XYZ
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr


__version__ = '1.0.0'

def main(dnaXYZ: str, monomerXYZ: str, atomGroups: str, out: str, distance: float, ignoreZeroPair: bool, **kwargs) -> None:

    with open(dnaXYZ) as dna_fh, open(monomerXYZ) as monomer_fh:
        allDistances = []
        while True:
            try:
                dna = read_XYZ(dna_fh)
                # Exclude type 3 to exclude 'inactive' TFs
                monomer = read_XYZ(monomer_fh, exclude='3')
                # Compute distance of monomers to each DNA bead
                result = cdist(dna['atoms'], monomer['atoms'], 'euclidean')
                # Find closest monomer to each bead
                minDistance = np.amin(result, axis=1) < distance
                allDistances.append(minDistance)
            except EOFError:
                break
        # Convert list of distance arrays to 2D numpy matrix
        allDistances = np.vstack(tuple(allDistances))

    # Read atomGroup and retrieve TU atom indexes
    atomGroupsDict = readJSON(atomGroups)
    TU_indexes = atomGroupsDict['TU']

    # Compute correlation ignore when both pairs are 0
    if ignoreZeroPair:
        correlations = []
        for idx1, col1 in enumerate(TU_indexes):
            for idx2, col2 in enumerate(TU_indexes):
                if col2 > col1: # Do not repeat duplicates
                    col1NoZero = allDistances[~((allDistances[:,idx1]==.0) & (allDistances[:,idx2]==.0)),idx1]
                    col2NoZero = allDistances[~((allDistances[:,idx1]==.0) & (allDistances[:,idx2]==.0)),idx2]
                    # Ensure col1 & col2 both have atleast 1 True & False (TF)
                    col1HaveTF = np.unique(col1NoZero).size == 2
                    col2HaveTF = np.unique(col2NoZero).size == 2
                    if col1HaveTF and col2HaveTF:
                        cor = np.corrcoef(col1NoZero, col2NoZero)[-1,0]
                        correlations.append([col1, col2, cor])
        correlations = pd.DataFrame(
            correlations, columns=['row', 'column', 'score'])
    else:
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
        pearson = pd.merge(correlations, pvalues)

    pearson.to_csv(out, index=False)


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
        'dnaXYZ', metavar='DNA',
        help='DNA coordinates in XYZ format')
    custom.add_argument(
        'monomerXYZ', metavar='MONOMER',
        help='Monomer coordinates in XYZ format')
    custom.add_argument(
        'atomGroups',
        help='Atom group assignments in JSON format.')
    custom.add_argument(
        '--distance', default=1.8, type=float,
        help='Distance threshold for active transcription (default: %(default)s)')
    custom.add_argument(
        '--out', default='TUcorrelation.csv',
        help='Contact matrix output (default: %(default)s)')
    custom.add_argument(
        '--ignoreZeroPair',  action='store_true',
        help='Ignore timepoint where both TU pairs are inactive when '
             'computing correlation. May be signicantly slower.')
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
