#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import numpy as np
from utilities import npz
import pyCommonTools as pct
from iced import normalization, filter
from timeit import default_timer as timer
from scipy.sparse import save_npz, load_npz, csc_matrix

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)
    parser.set_defaults(function=ICE)

    parser.add_argument(
        'matrix', type=npz, help='Input contact matrix')
    parser.add_argument(
        '-o', '--out', default='iced-contacts.npz', type=npz,
        help='Summed contact matrix (default: %(default)s)')
    parser.add_argument(
        '--filter_low', default=0.02, type=float,
        help='Percentage of lowest pairwise counts to filter.')

    return (pct.execute(parser))


def ICE(matrix: str, out: str, filter_low: float) -> None:

    log = pct.create_logger()
    start = timer()
    matrix = load_npz(matrix).toarray().astype(float)
    matrix = filter.filter_low_counts(matrix, percentage=filter_low)
    matrix = normalization.ICE_normalization(matrix)
    save_npz(out, csc_matrix(matrix))
    end = timer()
    log.info(f'Matrix ICED in {end - start} seconds.')


if __name__ == '__main__':
    sys.exit(main())
