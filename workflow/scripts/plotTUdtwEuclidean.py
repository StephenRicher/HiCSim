#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from typing import List
from utilities import coeff
import matplotlib.pyplot as plt


__version__ = '1.0.0'


def main(files: List, beadDistribution: str, out: str, minRep: int,
        fontSize: float, **kwargs) -> None:

    # Set global matplotlib fontisze
    plt.rcParams.update({'font.size': fontSize})

    allBeadDistribution = pd.read_csv(beadDistribution)
    polymerLength = len(allBeadDistribution.columns)

    TUdtw = pd.concat((pd.read_csv(file) for file in files))

    # Average across TU pairs
    TUdtw = TUdtw.groupby(
        ['TU1', 'TU2']).agg(['mean', 'count']).reset_index()

    # Set distance to np.nan for TU pairs with fewer than minRep samples
    TUdtw.loc[TUdtw[('distance', 'count')] < minRep, [('distance', 'mean')]] = np.nan

    # Convert to wide format
    TUdtw = TUdtw.pivot(index='TU1', columns='TU2', values=('distance', 'mean'))

    # Flip so diagonal is left to right
    TUdtw = TUdtw.iloc[::-1]

    fig, (ax1, ax2) = plt.subplots(
        2, gridspec_kw={'height_ratios': [6, 1]}, figsize=(8, 8))
    ax1 = sns.heatmap(TUdtw, square=True, vmin=0, cmap='viridis', ax=ax1)
    ax1.set_facecolor('xkcd:light grey')
    ax1.xaxis.set_label_text('')
    ax1.yaxis.set_label_text('')

    # Create 'n' xtick labels from start to end of polymer
    xticks = np.linspace(1, polymerLength, 10).astype(int)
    xticklabels = [i if i in xticks else '' for i in range(1, polymerLength + 1)]
    ax2 = sns.heatmap(allBeadDistribution,
        cmap='binary', vmin=0, vmax=len(files), ax=ax2,
        yticklabels=[''], xticklabels=xticklabels)

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
        '--out', default='TU-pairDistances.png',
        help='TU active pair distances heatmap (default: %(default)s)')
    custom.add_argument(
        '--fontSize', type=float, default=14,
        help='Font size for node name on circos plot (default: %(default)s)')
    custom.add_argument(
        '--minRep', type=int,
        help='Minimum replicates required for TU-TU interaction to be given '
             'included in the correlation/circos plots (default: %(default)s)')
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
