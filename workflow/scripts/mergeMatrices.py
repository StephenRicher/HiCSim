#!/usr/bin/env python3

""" Merge contact matrices """

import sys
import argparse
import numpy as np
from typing import List
from utilities import setDefaults, createMainParent
from scipy.sparse import save_npz, load_npz, csc_matrix


__version__ = '1.0.0'


def mergeMatrices(matrices: List, out: str, method: str):

    if method == 'median':
        average_matrix = compute_median(matrices)
    elif method == 'sum':
        average_matrix = compute_sum(matrices)
    else:
        average_matrix = compute_mean(matrices)

    save_npz(out, average_matrix)


def compute_median(matrices: List):
    """ Stack matrices along a third dimension and compute median """
    matrices = np.dstack([load_npz(matrix).toarray() for matrix in matrices])
    return csc_matrix(np.median(matrices, axis = 2))


def compute_sum(matrices: List):
    summed_matrix = 0
    for matrix in matrices:
        summed_matrix += load_npz(matrix)
    return summed_matrix


def compute_mean(matrices: List):
    summed_matrix = compute_sum(matrices)
    return summed_matrix / len(matrices)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=mergeMatrices)
    parser.add_argument('matrices', nargs='+', help='Input contact matrices')
    parser.add_argument(
        '--method', default='sum', choices=['mean', 'median', 'sum'],
        help='Method to compute average of matrices (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True, help='Summed contact matrix.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
