#!/usr/bin/env python3

""" Plot Radius of Gyration across replicates from LAMMPS output """


import sys
import argparse
import numpy as np
from scipy import stats
from typing import List
import matplotlib.pyplot as plt
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def plotRG(files: List, dpi: int, confidence: float, out: str):
    equilRG = files[0]
    postEquil = files[1:]

    data = np.dstack([np.loadtxt(rep) for rep in postEquil])
    mean = np.mean(data, axis=2)
    time = mean[:,0]
    mean = mean[:,1]

    # Filter duplicate timesteps from data
    time, uniqIdxs = np.unique(time, return_index=True)
    mean = mean[uniqIdxs]
    data = data[uniqIdxs,:,:]

    equilRG = np.loadtxt(equilRG)
    equilTime = equilRG[:,0]
    equilVal = equilRG[:,1]
    equilTime, uniqIdxs = np.unique(equilTime, return_index=True)
    equilVal = equilVal[uniqIdxs]

    allTime = np.concatenate([equilTime, time])
    allData = np.concatenate([equilVal, mean])

    n = len(postEquil)
    std_err = stats.sem(data, axis=2)[:,1]
    ci = std_err * stats.t.ppf((1 + confidence) / 2, n - 1)
    dydx = np.diff(allData)/ np.diff(allTime)
    midpoints = (allTime[1:] + allTime[:-1]) / 2

    fig, (ax1, ax2) = plt.subplots(2,1)

    ax1.plot(allTime, allData)
    ax1.fill_between(time, (mean - ci), (mean + ci), color='b', alpha=.1)
    ax1.axvline(x=1, color='Red', alpha=0.3)
    ax1.tick_params(labelsize=16)
    ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel('RG', fontsize=20)
    ax1.set_title(
        f'Radius of Gyration over time with {confidence:.0%} CI (n = {n})',
        loc='left', fontsize=20)

    ax2.plot(midpoints, dydx)
    ax2.axvline(x=1, color='Red', alpha=0.3)
    ax2.axhline(y=0, color='Red', alpha=0.3)
    ax2.tick_params(labelsize=16)
    ax2.set_xlabel('time', fontsize=20)
    ax2.set_ylabel('dRG/dt', fontsize=20)

    fig.set_size_inches(16, 14)
    fig.tight_layout()
    plt.savefig(fname=out, dpi=dpi)


def CI(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f'{x} not a floating-point literal.')
    if x <= 0.0 or x >= 1.0:
        raise argparse.ArgumentTypeError(f'{x} not in range (0.0, 1.0).')
    return x


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=plotRG)
    parser.add_argument(
        'files', nargs='*',
        help='Radius of Gyration output from LAMMPS.')
    parser.add_argument(
        '--dpi', type=int, default=600,
        help='Plot resolution (default: %(default)s)')
    parser.add_argument(
        '--confidence', type=CI, default=0.95,
        help='Confidence interval (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--out', required=True, help='Output plot filename.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
