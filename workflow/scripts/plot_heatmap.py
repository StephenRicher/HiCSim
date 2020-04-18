#!/usr/bin/env python3

""" Plot heatmap of contact frequency matrix """

import sys
import numpy as np
import seaborn as sns
import pyCommonTools as pct
import matplotlib.pyplot as plt


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)
    parser.set_defaults(function=plot_heatmap)

    parser.add_argument(
        'matrix', help='Input contact matrices')
    parser.add_argument(
        '--heatmap', default='heatmap.png',
        help='Output heatmap (default: %(default)s)')
    parser.add_argument(
        '--cmap', default='YlGn',
        help='Matplotlib colormap (default: %(default)s)')
    parser.add_argument(
        '--dpi', default=600, type=int,
        help='Heatmap resolution (default: %(default)s)')

    return (pct.execute(parser))


def plot_heatmap(matrix: str, heatmap: str, dpi: int, cmap: str) -> None:

    matrix = np.loadtxt(matrix)
    sns.heatmap(np.log10(matrix + 1), cmap)
    plt.savefig(heatmap, dpi=dpi)


if __name__ == '__main__':
    sys.exit(main())
