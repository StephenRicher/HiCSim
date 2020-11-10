#!/usr/bin/env python3

""" Plot and cluster TU activation over time. """

import sys
import logging
import argparse
import pandas as pd
import seaborn as sns
from typing import List
import matplotlib.pyplot as plt
from utilities import setDefaults

__version__ = '1.0.0'


def plotTUactivation(files: List, out: str):

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


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    custom.add_argument(
        'files', nargs='+', help='Input TF activation tables.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True, help='Output plot filename.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(plotTUactivation(**vars(args)))
