#!/usr/bin/env python3

""" Plot heatmap of contact frequency matrix """

import sys
import numpy as np
import seaborn as sns
import pyCommonTools as pct
import matplotlib.pyplot as plt
from scipy.sparse import load_npz


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)
    parser.set_defaults(function=plot_heatmap)

    parser.add_argument(
        'matrix', help='Input contact matrices')
    parser.add_argument(
        '--transform', default='none', choices=['none', 'log10', 'obsexp'],
        help='Transformation to apply to counts (default: %(default)s)')
    parser.add_argument(
        '--heatmap', default='heatmap.png',
        help='Output heatmap (default: %(default)s)')
    parser.add_argument(
        '--cmap', default='YlGn',
        help='Matplotlib colormap (default: %(default)s)')
    parser.add_argument(
        '--vmin', default=None, type=float,
        help='Min value to anchor colourmap (default: %(default)s)')
    parser.add_argument(
        '--vmax', default=None, type=float,
        help='Max value to anchor colourmap (default: %(default)s)')
    parser.add_argument(
        '--dpi', default=600, type=int,
        help='Heatmap resolution (default: %(default)s)')

    return (pct.execute(parser))


def plot_heatmap(
    matrix: str, transform: str, heatmap: str, dpi: int,
    cmap: str, vmin: float, vmax: float) -> None:

    matrix = load_npz(matrix).todense()

    if transform == 'log10':
        matrix = np.log10(matrix + 1)
    elif transform == 'obsexp':
        matrix = obsexp(matrix)

    sns.heatmap(matrix, cmap=cmap, vmax=2)
    plt.savefig(heatmap, dpi=dpi)


def obsexp(matrix):
    nbins = len(matrix)
    for offset in range(-nbins + 1, nbins):
        len_diag = nbins - abs(offset)
        sum_diag = np.trace(matrix, offset = offset)
        expected = sum_diag / len_diag
        if expected == 0:
            matrix[kth_diag_indices(matrix, offset)] = np.nan
        else:
            matrix[kth_diag_indices(matrix, offset)] /= expected
    return matrix


def kth_diag_indices(a, k):
    """ Retrieve indicies of diagonal at offset k """

    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


if __name__ == '__main__':
    sys.exit(main())
