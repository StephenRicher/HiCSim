#!/usr/bin/env python3

""" Calculate TU-TU correlation """

import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from utilities import pearsonr_pval, countPair
from statsmodels.stats.multitest import fdrcorrection

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
    allValues['p(adj)'] = fdrcorrection(allValues['p'])[1]
    allValues.to_csv(out)


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
