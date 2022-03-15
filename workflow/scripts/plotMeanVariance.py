#!/usr/bin/env python3

""" Plot mean-variance relationship of TU activity expression """

import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def plotMeanVariance(infile: str, out: str, fontSize: float) -> None:

    # Set global matplotlib fontisze
    plt.rcParams.update({'font.size': fontSize})
    stats = pd.read_csv(infile, usecols=['TU', 'activity'])
    stats = stats.groupby('TU').agg(['mean', 'var'])

    fig, ax = plt.subplots()
    stats.plot.scatter(
        x=('activity', 'mean'), y=('activity',  'var'),
        c=stats.index, ax=ax, cmap='viridis')

    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches='tight')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=plotMeanVariance)
    parser.add_argument(
        'infile', default=sys.stdin, help='Input TU stats (default: stdin)')
    parser.add_argument(
        '--fontSize', type=float, default=14,
        help='Font size for plot (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True, help='Output plot filename.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
