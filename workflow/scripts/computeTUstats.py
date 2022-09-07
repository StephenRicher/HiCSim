#!/usr/bin/env python3

""" Compute TU expression metrics """


import sys
import zlib
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict
from utilities import readJSON, pearsonr_pval, countPair


def main(infile: str, out: str, TUpairStats: str) -> None:

    TUinfo = pd.read_csv(infile, usecols=['timestep', 'ID', 'active'])
    TUinfo = TUinfo.sort_values('timestep').groupby('ID')
    TUstats = defaultdict(list)
    TUactivity = {}

    for TU, info in TUinfo:
        activeBursts, inactiveBursts = getConsecutiveBool(info['active'])
        totalActivity = sum(info['active'])
        # Calculate timestep interval
        simTime = info['timestep'].max()
        nSteps = len(info['active'])
        timestep = simTime / nSteps

        # Compute minimum compression ratio for given length
        minCompress = len(zlib.compress(b'.' * nSteps)) / len(info['active'])
        TUstats['TU'].append(TU)
        TUstats['simTime'].append(simTime)
        TUstats['nSteps'].append(nSteps)
        # Number of consecutive ON periods per unit time 10e6
        TUstats['nBursts'].append(len(activeBursts))
        TUstats['nBursts / 10e6'].append((len(activeBursts) / simTime) * 10e6)
        TUstats['activity'].append(totalActivity / nSteps)
        TUstats['burstLengthMean'].append(np.mean(activeBursts) * timestep)
        TUstats['burstLengthStd'].append(np.std(activeBursts) * timestep)
        TUstats['inactiveLengthMean'].append(np.mean(inactiveBursts) * timestep)
        TUstats['inactiveLengthStd'].append(np.std(inactiveBursts) * timestep)
        # Convert active boolean array to binary string
        boolString = ''.join(str(int(x)) for x in np.array(info['active']))
        activeCompressed = zlib.compress(boolString.encode("utf-8"))
        TUstats['complexity'].append(
            (len(activeCompressed) / len(info['active'])) / minCompress)
        # Store per-TU activity across all timepoints in seperate dict
        TUactivity[TU] = info['active'].to_list()

    TUstats = pd.DataFrame(TUstats)
    TUstats.to_csv(out, index=False)

    TUactivity = pd.DataFrame(TUactivity)
    TUcorr = TUactivity.corr('pearson').stack().rename('r')
    TUpval = TUactivity.corr(method=pearsonr_pval).stack().rename('p')
    TUactivity = (
        pd.merge(TUcorr, TUpval, left_index=True, right_index=True)
        .reset_index()
        .rename({'level_0': 'row', 'level_1': 'column'}, axis=1)
    )
    TUactivity.to_csv(TUpairStats, index=False)


def getConsecutiveBool(arr):
    """ Return lengths of Boolean chunks in boolean array """

    arr = np.array(arr)
    # Compute boundaries
    boundaries = np.concatenate(([arr[0]], arr[:-1] != arr[1:], [True]))
    # Retrieve index positions of boundaries
    boundaries  = np.where(boundaries)[0]
    # Get differences between indexes
    differences = np.diff(boundaries)
    truthChunks = differences[::2]
    falseChunks = differences[1::2]

    # Edge case where array start with False, first chunk is missed.
    if arr[0] == False:
        # First index of boundaries corresponds to length of first Falses
        falseChunks = np.insert(falseChunks, 0, boundaries[0])
    return truthChunks, falseChunks


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infile', default=sys.stdin,
        help='Input TU activation table.')
    parser.add_argument('--out', default=sys.stdout,
        help='File to save TU expression metrics.')
    parser.add_argument('--TUpairStats', default=sys.stderr,
        help='Optionally output TU-TU correlation with '
             'relative positioning info.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
