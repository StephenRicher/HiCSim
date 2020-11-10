#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import logging
import argparse
import numpy as np
from typing import List
from utilities import npz
from scipy.sparse import save_npz, load_npz, csc_matrix


__version__ = '1.0.0'


def main(matrices: List, out: str, method: str, **kwargs) -> None:

    if method == 'median':
        average_matrix = compute_median(matrices)
    elif method == 'sum':
        average_matrix = compute_sum(matrices)
    else:
        average_matrix = compute_mean(matrices)

    save_npz(out, average_matrix)


def compute_median(matrices: List):
    # Stack matrices along a third dimension and compute median
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


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'matrices', nargs='+',
        help='Input contact matrices')
    custom.add_argument(
        '--out', default='averaged-contacts.npz', type=npz,
        help='Summed contact matrix (default: %(default)s)')
    custom.add_argument(
        '--method', default='sum', choices=['mean', 'median', 'sum'],
        help='Method to compute average of matrices (default: %(default)s)')
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
