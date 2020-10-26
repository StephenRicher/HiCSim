#!/usr/bin/env python3

""" Plot mean-variance relationship of TU activity expression """

import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(infile: str, out: str, fontSize: float) -> None:

    # Set global matplotlib fontisze
    plt.rcParams.update({'font.size': fontSize})
    stats = pd.read_csv(infile, usecols=['TU','activity'])
    stats = stats.groupby('TU').agg(['mean', 'var'])

    fig, ax = plt.subplots()
    stats.plot.scatter(x=('activity', 'mean'), y=('activity',  'var'),
        c=stats.index, ax=ax, cmap='viridis')

    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches='tight')


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infile', default=sys.stdin,
        help='Input TU stats.')
    parser.add_argument('out',
        help='File to save scatter plot.')
    parser.add_argument(
        '--fontSize', type=float, default=14,
        help='Font size for node name on circos plot (default: %(default)s)')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
