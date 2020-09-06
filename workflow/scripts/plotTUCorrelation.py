#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import logging
import argparse
import pandas as pd
import seaborn as sns
from typing import List
from matplotlib import pyplot as plt

__version__ = '1.0.0'


def main(files: List, out: str, **kwargs) -> None:

    # Read all per-replicates correlation dataframes into 1 dataframe
    correlation = pd.concat((pd.read_csv(file) for file in files))

    # Get pairwise mean correlation across replicates
    correlation = correlation.groupby(['row', 'column']).mean()

    # Convert to wide format matrix
    correlation = correlation.reset_index().pivot(
        index='row', columns='column', values='score')

    # Flip vertically to ensure diagonal goes from bottom left to top right
    correlation = correlation.iloc[::-1]
    #vmin=-0.3, vmax=0.3
    ax = sns.heatmap(correlation, cmap='bwr', center=0)
    plt.savefig(out)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'files', nargs='+',
        help='Input TF distance tables.')
    custom.add_argument(
        '--out', default='TF-correlation.png',
        help='TF contact correlation heatmap (default: %(default)s)')
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
