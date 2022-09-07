#!/usr/bin/env python3

""" Plot pairwise TU correlation and heatmap """

import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def main(file: str, beadDistribution: str, meanHeatmap: str, sumHeatmap: str,
        circos: str, minRep: int, nReps: int, pvalue: float,
        vmin: float, vmax: float, fontSize: float, **kwargs) -> None:

    # Set global matplotlib fontisze
    plt.rcParams.update({'font.size': fontSize})

    correlation = pd.read_csv(file)

    allBeadDistribution = pd.read_csv(beadDistribution)
    polymerLength = len(allBeadDistribution.columns)

    nodeNames = correlation['row'].unique()

    for method in ['r', 'count']:
        transform = correlation
        if method == 'r':
            transform = transform.set_index(['row', 'column'])
            # Set values to NaN if replicate count less than threshold
            transform.loc[transform['count'] < minRep, ['r', 'p(adj)']] = np.nan

            G = nx.Graph()
            G.add_nodes_from(sorted(nodeNames))
            colours = []
            for row, column in transform.index:
                if row >= column:
                    continue
                info = transform.loc[row, column]
                if info['p(adj)'] >= pvalue:
                    continue
                colours.append(info['r'])
                G.add_edge(row, column)
            nx.draw_circular(G,node_color='White',
                         node_size=0,
                         font_size=fontSize,
                         edge_color=colours,
                         edge_cmap=plt.cm.RdBu_r,
                         edge_vmin=vmin, edge_vmax=vmax,
                         with_labels=True)
            plt.savefig(circos, bbox_inches='tight')

            transform = transform.reset_index()

        # Set diagonal to NaN to hide trivial auto-correlation
        transform.loc[transform.row == transform.column, method] = np.nan
        transform = transform.pivot(index='row', columns='column', values=method)

        # Flip vertically to ensure diagonal goes from bottom left to top right
        transform = transform.iloc[::-1]
        fig, (ax1, ax2) = plt.subplots(2, gridspec_kw={'height_ratios': [6, 1]}, figsize=(8, 8))

        if method == 'r':
            ax1 = sns.heatmap(transform, square=True,
                cmap='bwr', center=0, vmin=vmin, vmax=vmax, ax=ax1)
            # Ensure masked cells are not within 'bwr' colour map.
            ax1.set_facecolor('xkcd:light grey')
            out = meanHeatmap
        else:
            ax1 = sns.heatmap(transform, square=True,
                cmap='binary', vmin=0, vmax=nReps, ax=ax1)
            out = sumHeatmap

        ax1.xaxis.set_label_text('')
        ax1.yaxis.set_label_text('')

        # Create 'n' xtick labels from start to end of polymer
        xticks = np.linspace(1, polymerLength, 10).astype(int)
        xticklabels = [i if i in xticks else '' for i in range(1, polymerLength + 1)]
        ax2 = sns.heatmap(allBeadDistribution,
            cmap='binary', vmin=0, vmax=nReps, ax=ax2,
            yticklabels=[''], xticklabels=xticklabels)

        # Unset xtick markers for empty string xlabels
        for index, xtick in enumerate(ax2.xaxis.get_major_ticks()):
            if xticklabels[index] == '':
                xtick.set_visible(False)

        fig.tight_layout()
        fig.savefig(out, dpi=300, bbox_inches='tight')


def coeff(value):
    ivalue = float(value)
    if not -1 <= ivalue <= 1:
        raise argparse.ArgumentTypeError(f'{value} must be between -1 and 1.')
    return ivalue


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=main)
    parser.add_argument('beadDistribution',
        help='Bead distribution output of writeTUdistribution.')
    parser.add_argument(
        'file',
        help='Input correlation table.')
    parser.add_argument(
        '--meanHeatmap', default='TF-meanCorrelation.png',
        help='TF mean contact correlation heatmap (default: %(default)s)')
    parser.add_argument(
        '--sumHeatmap', default='TF-replicateCount.png',
        help='Replicate count heatmap for each pairwise correlation (default: %(default)s)')
    parser.add_argument(
        '--circos', default='TU-CircosPlot.png',
        help='TF circos plot for signicant correlations (default: %(default)s)')
    parser.add_argument(
        '--fontSize', type=float, default=14,
        help='Font size for node name on circos plot (default: %(default)s)')
    parser.add_argument(
        '--pvalue', type=float, default=10**-6,
        help='P-value threshold for filtering TU correlations before '
             'plotting on circos plot (default: %(default)s)')
    parser.add_argument(
        '--minRep', type=int, default=1,
        help='Minimum replicates required for TU-TU interaction to be given '
             'included in the correlation/circos plots (default: %(default)s)')
    parser.add_argument(
        '--nReps', type=int,
        help='Number of reps in dataset - used to scale colour scheme of '
             'count matrix')
    parser.add_argument(
        '--vmin', type=coeff, default=-0.3,
        help='Minimum value of colour scale. (default: %(default)s)')
    parser.add_argument(
        '--vmax', type=coeff, default=0.3,
        help='Maximum value of colour scale. (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
