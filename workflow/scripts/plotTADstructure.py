#!/usr/bin/env python3

""" Plot TAD structure over time """


import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utilities import setDefaults


__version__ = '1.0.0'


def plotTADstructure(beadTADstatus: str, dpi: int, out: str):

    dat = pd.read_csv(beadTADstatus)
    data = dat.pivot(columns='time', index='ID', values='TADstatus')

    fig, ax = plt.subplots()
    statuses = dat['TADstatus'].unique()
    cmap = sns.mpl_palette('Dark2', len(statuses) - 2)
    cmap.insert(0, ('#d3d3d3'))
    cmap.insert(0, ('#000000'))
    sns.heatmap(data, cmap=cmap, xticklabels='auto', cbar=False,
                vmin=min(statuses)-0.5, vmax=max(statuses)+0.5, ax=ax)
    ax.set_xlabel('Time')
    ax.set_ylabel('Bead ID')
    fig.set_size_inches(16, 9)
    fig.tight_layout()
    fig.savefig(fname=out, dpi=dpi)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'beadTADstatus', help='TADstatus output of runLammpsSimulation.')
    parser.add_argument(
        '--dpi', type=int, default=600,
        help='Plot resolution (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True, help='Output plot filename.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(plotTADstructure(**vars(args)))
