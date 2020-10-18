#!/usr/bin/env python3

""" Plot and cluster TU activation over time. """

import sys
import logging
import argparse
import pandas as pd
import seaborn as sns
from typing import List
import matplotlib.pyplot as plt

__version__ = '1.0.0'


def main(files: List, out: str, **kwargs) -> None:

    activations = pd.concat((pd.read_csv(file) for file in files))
    activations = activations.groupby(['time', 'id']).sum().reset_index()
    activations = activations.pivot(
        index='id', columns='time', values='active')

    #plot = sns.clustermap(activations, cmap='binary', col_cluster=False, vmin=0, vmax=len(files))
    #plot.savefig(out, dpi=300, bbox_inches='tight')

    fig, ax = plt.subplots()
    ax = sns.heatmap(activations, cmap='binary', vmin=0, vmax=len(files), ax=ax)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches='tight')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'files', nargs='+',
        help='Input TF activation tables.')
    custom.add_argument(
        '--out', default='TU-activation.png',
        help='TF mean contact correlation heatmap (default: %(default)s)')

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
