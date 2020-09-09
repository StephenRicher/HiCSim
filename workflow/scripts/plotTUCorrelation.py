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
import matplotlib.pyplot as plt
from nxviz.plots import CircosPlot
#from matplotlib import pyplot as plt

__version__ = '1.0.0'


def main(files: List, heatmap: str, circos: str, pvalue: float, **kwargs) -> None:

    # Read all per-replicates correlation dataframes into 1 dataframe
    correlation = pd.concat((pd.read_csv(file) for file in files))

    nodeNames = correlation['row'].unique()

    # Get pairwise mean correlation across replicates
    correlation = correlation.groupby(['row', 'column']).mean()

    filtered = correlation[correlation['p-value'] < pvalue].reset_index()
    filtered['direction'] = np.where(filtered['r'] > 0, 'Up', 'Down')
    filtered['width'] = 3
    G = nx.from_pandas_edgelist(
        filtered.reset_index(), "row", "column", ['direction', 'width'])
    nx.set_node_attributes(G, pd.Series(nodeNames, index=nodeNames), 'name')
    c = CircosPlot(G,
               edge_color="direction",
               edge_width="width",
               node_labels=True,
               node_order='name',
               edgeprops = {"facecolor": "none", "alpha": 1})
    c.draw()
    plt.tight_layout()
    plt.savefig(circos)

    correlation = correlation.reset_index()
    # Set diagonal to 0 to hide trivial auto-correlation
    correlation.loc[correlation.row == correlation.column, 'r'] = 0
    correlation = correlation.pivot(index='row', columns='column', values='r')

    # Flip vertically to ensure diagonal goes from bottom left to top right
    correlation = correlation.iloc[::-1]
    fig, ax = plt.subplots()
    ax = sns.heatmap(correlation, cmap='bwr', center=0, vmin=-0.3, vmax=0.3)
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
