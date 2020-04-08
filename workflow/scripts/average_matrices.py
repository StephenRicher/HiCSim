#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import numpy as np
from typing import List
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)
    parser.set_defaults(function=average_heatmap)

    parser.add_argument(
        'matrices', nargs='+',
        help='Input contact matrices')
    parser.add_argument(
        '-o', '--out', default='contacts-summed.txt',
        help='Summed contact matrix (default: %(default)s)')
    parser.add_argument(
        '--method', default='median', choices=['mean', 'median'],
        help='Method to compute average of matrices (default: %(default)s)')

    return (pct.execute(parser))


def compute_median(matrices: List):
    # Stack matrices along a third dimension and compute median
    matrices = np.dstack([np.loadtxt(matrix) for matrix in matrices])

    return np.median(matrices, axis = 2)


def compute_mean(matrices: List):
    summed_matrix = 0
    for matrix in matrices:
        summed_matrix += np.loadtxt(matrix)

    return summed_matrix / len(matrices)


def average_heatmap(matrices: List, out: str, method: str) -> None:

    if method == 'median':
        average_matrix = compute_median(matrices)
    else:
        average_matrix = compute_mean(matrices)

    np.savetxt(out, average_matrix)


if __name__ == '__main__':
    sys.exit(main())
