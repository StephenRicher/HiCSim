#!/usr/bin/env python3

""" Compute TU expression metrics """


import sys
import zlib
import logging
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict


__version__ = '1.0.0'


def main(file: str, outdir: str, **kwargs) -> None:

    TUinfo = pd.read_csv(file, usecols=['time', 'id', 'active'])
    TUinfo = TUinfo.sort_values('time').groupby('id')

    TUstats = defaultdict(list)

    for TU, info in TUinfo:
        activeBursts, inactiveBursts = getConsecutiveBool(info.active)
        totalActivity = sum(activeBursts)

        # Calculate timestep interval
        simTime = info.time.max()
        nSteps = len(info.active)
        timestep = simTime / nSteps

        # Compute minimum compression ratio for given length
        minCompress = len(zlib.compress(b'.' * nSteps)) / len(info.active)

        TUstats['TU'].append(TU)
        TUstats['simTime'].append(simTime)
        TUstats['nSteps'].append(nSteps)
        # Number of consecutive ON periods per unit time 10e6
        TUstats['nBursts'].append(len(activeBursts))
        TUstats['nBursts / 10e6'].append((len(activeBursts) / simTime) * 10e6)
        TUstats['activity'].append(totalActivity / timepoints)
        TUstats['burstLengthMean'].append(np.mean(activeBursts) * timestep)
        TUstats['burstLengthStd'].append(np.std(activeBursts) * timestep)
        TUstats['inactiveLengthMean'].append(np.mean(inactiveBursts) * timestep)
        TUstats['inactiveLengthStd'].append(np.std(inactiveBursts) * timestep)

        # Convert active boolean array to binary string
        boolString = ''.join(str(int(x)) for x in np.array(info.active))
        activeCompressed = zlib.compress(boolString.encode("utf-8"))
        TUstats['complexity'].append(
            (len(activeCompressed) / len(info.active)) / minCompress)


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

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('file', default=sys.stdin,
        help='Input TU activation table.')
    custom.add_argument('--out', default=sys.stdout,
        help='FIle to save TU expression metrics.')


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
