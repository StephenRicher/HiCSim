#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import numpy as np
from typing import List
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)
    parser.set_defaults(function=plot_heatmap)

    parser.add_argument(
        'matrices', nargs='+',
        help='Input contact matrices')
    parser.add_argument(
        '-o', '--out', default='contacts-summed.txt',
        help='Summed contact matrix (default: %(default)s)')

    return (pct.execute(parser))


def plot_heatmap(matrices: List, out: str) -> None:

    summed_matrix = 0
    for matrix in matrices:
        summed_matrix += np.loadtxt(matrix)
    average_matrix = summed_matrix / len(matrices)

    np.savetxt(out, average_matrix)


if __name__ == '__main__':
    sys.exit(main())
