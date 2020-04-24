#!/usr/bin/env python3

""" Plot heatmap of contact frequency matrix """

import sys
import numpy as np
import seaborn as sns
from typing import Dict
from bisect import bisect
import pyCommonTools as pct
from utilities import region, scale
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
        '--vmin', default=None, type=scale,
        help='Min value to anchor colourmap (default: %(default)s)')
    parser.add_argument(
        '--vmax', default=None, type=scale,
        help='Max value to anchor colourmap (default: %(default)s)')
    parser.add_argument(
        '--dpi', default=600, type=int,
        help='Heatmap resolution (default: %(default)s)')
    parser.add_argument(
        '--region', type=region,
        metavar='START-END',
        help='Specify bin size, start and end coordinates of contact matrix '
        'for plotting axis labels.')
    parser.add_argument(
        '--viewpoint', type=region,
        metavar='START-END',
        help='Specify start and end viewpoint of contact matrix. '
        'Must also specify --region argument.')

    return (pct.execute(parser))


def get_region_coordinates(region, size):
    coordinates = np.linspace(region['start'], region['end'], size)
    return coordinates.astype(int)


def trim_duplicates(values):
    """ Remove duplicates to avoid duplicate labels in heatmap """

    deduplicated = []
    for i, value in enumerate(values):
        # Exclude first value as coordinate probably not exact
        if i == 0:
            first = value
        if value not in deduplicated and value != first:
            deduplicated.append(value)
        else:
            deduplicated.append('')
    return deduplicated


def validate_viewpoints(region, viewpoint):
    """ Ensure viewpoint coordinates are within region coordinates. """

    log = pct.create_logger()
    if (viewpoint['start'] < region['start']
            or viewpoint['start'] > region['end']):
        log.error(
            f'Start viewpoint {viewpoint["start"]} is not within '
            f'region {region["start"]}-{region["end"]}.')
        sys.exit(1)
    if (viewpoint['end'] < region['start']
            or viewpoint['end'] > region['end']):
        log.error(
            f'End viewpoint {viewpoint["end"]} is not within '
            f'region {region["start"]}-{region["end"]}.')
        sys.exit(1)


def plot_heatmap(
    matrix: str, transform: str, heatmap: str, dpi: int, cmap: str,
    vmin: float, vmax: float, viewpoint: Dict, region: Dict) -> None:

    matrix = load_npz(matrix).todense()
    if region:
        region_bins = get_region_coordinates(region, len(matrix))
        if viewpoint:
            validate_viewpoints(region, viewpoint)
            start_view = bisect(region_bins, viewpoint['start'])
            end_view = bisect(region_bins, viewpoint['end'])
            region_bins = region_bins[start_view:end_view]
            matrix = matrix[start_view:end_view,start_view:end_view]
        if region_bins[-1] - region_bins[0] > 100000:
            divisor = 1000000
            label = 'Mb'
        else:
            divisor = 100000
            label = 'Kb'
        region_bins = [f'{bin/divisor:.1f}' for bin in region_bins]
        region_bins = trim_duplicates(region_bins)
    else:
        region_bins = None

    if transform == 'log10':
        matrix = np.log10(matrix + 1)
    elif transform == 'obsexp':
        matrix = obsexp(matrix)

    ax = sns.heatmap(matrix, cmap=cmap, vmin=vmin, vmax=vmax,
        xticklabels=region_bins, yticklabels=region_bins)
    remove_empty_ticks(ax)

    plt.yticks(rotation=0, fontsize='small')
    plt.xticks(rotation=45, fontsize='small')
    plt.xlabel(label)
    plt.ylabel(label)
    plt.tight_layout()
    plt.savefig(heatmap, dpi=dpi)


def remove_empty_ticks(ax):
    """ Remove ticks from plot that are labelled as empty strings """

    xticks = ax.xaxis.get_major_ticks()
    for i, label in enumerate(ax.get_xticklabels()):
        if label.get_text() == '':
            xticks[i].set_visible(False)

    yticks = ax.yaxis.get_major_ticks()
    for i, label in enumerate(ax.get_yticklabels()):
        if label.get_text() == '':
            yticks[i].set_visible(False)


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
