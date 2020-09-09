#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
from typing import List
from utilities import coeff
import matplotlib.pyplot as plt


__version__ = '1.0.0'


def main(files: List, heatmap: str, circos: str, pvalue: float, vmin: float, vmax: float, **kwargs) -> None:

    # Read all per-replicates correlation dataframes into 1 dataframe
    correlation = pd.concat((pd.read_csv(file) for file in files))

    nodeNames = correlation['row'].unique()

    # Get pairwise mean correlation across replicates
    correlation = correlation.groupby(['row', 'column']).mean()

    G = nx.Graph()
    G.add_nodes_from(sorted(nodeNames))
    colours = []
    for row, column in correlation.index:
        if row >= column:
            continue
        info = correlation.loc[row, column]
        if info['p-value'] >= pvalue:
            continue
        colours.append(info['r'])
        G.add_edge(row, column)
    nx.draw_circular(G,node_color='White',
                 node_size=0,
                 font_size=14,
                 edge_color=colours,
                 edge_cmap=plt.cm.RdBu_r,
                 edge_vmin=vmin, edge_vmax=vmax,
                 with_labels=True)
    plt.tight_layout()
    plt.savefig(circos)

    correlation = correlation.reset_index()
    # Set diagonal to 0 to hide trivial auto-correlation
    correlation.loc[correlation.row == correlation.column, 'r'] = 0
    correlation = correlation.pivot(index='row', columns='column', values='r')

    # Flip vertically to ensure diagonal goes from bottom left to top right
    correlation = correlation.iloc[::-1]
    fig, ax = plt.subplots()
    ax = sns.heatmap(correlation, cmap='bwr', center=0, vmin=vmin, vmax=vmax)
    fig.tight_layout()
    fig.savefig(heatmap)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'files', nargs='+',
        help='Input TF distance tables.')
    custom.add_argument(
        '--heatmap', default='TF-Correlation.png',
        help='TF contact correlation heatmap (default: %(default)s)')
    custom.add_argument(
        '--circos', default='TU-CircosPlot.png',
        help='TF circos plot for signicant correlations (default: %(default)s)')#
    custom.add_argument(
        '--pvalue', type=float, default=10**-6,
        help='P-value threshold for filtering TU correlations before '
             'plotting on circos plot (default: %(default)s)')
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
