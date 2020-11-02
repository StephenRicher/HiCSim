#!/usr/bin/env python3

""" Calculate TU-TU correlation """

import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


__version__ = '1.0.0'

def main(infile: str, out: str) -> None:

    activityMatrix = pd.read_csv(infile).pivot(
        index='rep', columns='TU', values='activity')

    allValues = []
    for method in ['pearson', pearsonr_pval, countPair]:
        values = activityMatrix.corr(method=method).stack()
        allValues.append(values)
    allValues = pd.concat(allValues, axis=1)
    allValues.index.names = ['row', 'column']
    allValues.columns = ['r', 'p','count']

    allValues.to_csv(out)


def pearsonr_pval(x,y):
    try:
        return pearsonr(x,y)[1]
    except ValueError:
        return np.nan


def countPair(x, y):
    """ Return count of valid pairs (both not nan) """

    # Indices where both x and y are NOT np.nan
    validIndices = np.intersect1d(
        np.where(~np.isnan(x)),
        np.where(~np.isnan(y)))
    return len(validIndices)


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infile', default=sys.stdin,
        help='Input TU stats.')
    parser.add_argument('--out', default=sys.stdout,
        help='File to save correlation table.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
