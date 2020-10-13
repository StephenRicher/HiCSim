#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import logging
import argparse
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
from typing import List
from utilities import coeff
import matplotlib.pyplot as plt


__version__ = '1.0.0'


def main(files: List, beadDistribution: str, meanHeatmap: str, sumHeatmap: str,
        circos: str, minRep: int, pvalue: float, vmin: float, vmax: float,
        fontSize: float, **kwargs) -> None:

    # Read all per-replicates correlation dataframes into 1 dataframe
    correlation = pd.concat((pd.read_csv(file) for file in files))

    allBeadDistribution = pd.read_csv(beadDistribution)
    polymerLength = len(allBeadDistribution.columns)

    nodeNames = correlation['row'].unique()

    # Get pairwise mean correlation across replicates
    correlation = correlation.groupby(['row', 'column']).agg(['mean', 'count'])

    for method in ['mean', 'count']:
        transform = correlation
        if method == 'mean':
            # Set mean values to NaN if replicate count less than threshold
            if minRep:
                transform.loc[transform[('r','count')] < minRep, [('r', 'mean'), ('p-value','mean')]] = np.nan

            G = nx.Graph()
            G.add_nodes_from(sorted(nodeNames))
            colours = []
            for row, column in transform.index:
                if row >= column:
                    continue
                info = transform.loc[row, column]
                if info['p-value']['mean'] >= pvalue:
                    continue
                colours.append(info['r']['mean'])
                G.add_edge(row, column)
            nx.draw_circular(G,node_color='White',
                         node_size=0,
                         font_size=fontSize,
                         edge_color=colours,
                         edge_cmap=plt.cm.RdBu_r,
                         edge_vmin=vmin, edge_vmax=vmax,
                         with_labels=True)
            plt.savefig(circos, bbox_inches='tight')

        # Extract the 'r' and 'p-value' columns corresponding to method
        transform = transform.iloc[:, transform.columns.get_level_values(1)==method]
        # Remove redundant 'method' multi-level index
        transform.columns = transform.columns.droplevel(1)

        transform = transform.reset_index()

        # Set diagonal to 0 to hide trivial auto-correlation
        transform.loc[transform.row == transform.column, 'r'] = np.nan
        transform = transform.pivot(index='row', columns='column', values='r')

        # Flip vertically to ensure diagonal goes from bottom left to top right
        transform = transform.iloc[::-1]
        fig, (ax1, ax2) = plt.subplots(2, gridspec_kw={'height_ratios': [6, 1]}, figsize=(8, 8))

        if method == 'mean':
            ax1 = sns.heatmap(transform, square=True,
                cmap='bwr', center=0, vmin=vmin, vmax=vmax, ax=ax1)
            # Ensure masked cells are not within 'bwr' colour map.
            ax1.set_facecolor('xkcd:light grey')
            out = meanHeatmap
        else:
            ax1 = sns.heatmap(transform, square=True,
                cmap='binary', vmin=0, vmax=len(files), ax=ax1)
            out = sumHeatmap

        ax1.xaxis.set_label_text('')
        ax1.yaxis.set_label_text('')

        # Create 'n' xtick labels from start to end of polymer
        xticks = np.linspace(1, polymerLength, 10).astype(int)
        xticklabels = [i if i in xticks else '' for i in range(1, polymerLength + 1)]
        ax2 = sns.heatmap(allBeadDistribution,
            cmap='binary', vmin=0, vmax=len(files), ax=ax2,
            yticklabels=[''], xticklabels=xticklabels)

        # Unset xtick markers for empty string xlabels
        for index, xtick in enumerate(ax2.xaxis.get_major_ticks()):
            if xticklabels[index] == '':
                xtick.set_visible(False)

        fig.tight_layout()
        fig.savefig(out, dpi=300, bbox_inches='tight')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'beadDistribution', help='Bead distribution output of writeTUdistribution.')
    custom.add_argument(
        'files', nargs='+', help='Input TF distance tables.')
    custom.add_argument(
        '--meanHeatmap', default='TF-meanCorrelation.png',
        help='TF mean contact correlation heatmap (default: %(default)s)')
    custom.add_argument(
        '--sumHeatmap', default='TF-replicateCount.png',
        help='Replicate count heatmap for each pairwise correlation (default: %(default)s)')
    custom.add_argument(
        '--circos', default='TU-CircosPlot.png',
        help='TF circos plot for signicant correlations (default: %(default)s)')
    custom.add_argument(
        '--fontSize', type=float, default=14,
        help='Font size for node name on circos plot (default: %(default)s)')
    custom.add_argument(
        '--pvalue', type=float, default=10**-6,
        help='P-value threshold for filtering TU correlations before '
             'plotting on circos plot (default: %(default)s)')
    custom.add_argument(
        '--minRep', type=int,
        help='Minimum replicates required for TU-TU interaction to be given '
             'included in the correlation/circos plots (default: %(default)s)')
    custom.add_argument(
        '--vmin', type=coeff, default=-0.3,
        help='Minimum value of colour scale. (default: %(default)s)')
    custom.add_argument(
        '--vmax', type=coeff, default=0.3,
        help='Maximum value of colour scale. (default: %(default)s)')
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
