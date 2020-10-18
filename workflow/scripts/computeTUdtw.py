#!/usr/bin/env python3

""" Calculate Euclidean distance between aligned time series using DTW. """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
from dtaidistance import dtw_ndim

__version__ = '1.0.0'

def main(file: str, out: str, downsample: int, **kwargs) -> None:

    # Read input and sort by TU then by time and group by TU
    info = pd.read_csv(file)

    if downsample:
        info = info[info.time % downsample == 0]

    names = info['id'].unique()
    info = info.sort_values(['id', 'time']).groupby(['id'])

    # Loop through TUs and add to timeseries data structure
    allPositions = []
    for name, timeseries in info:
        allPositions.append(timeseries[['xPos', 'yPos', 'zPos']].to_numpy())

    # Compute time aligned Euclidean distance
    dtw = dtw_ndim.distance_matrix_fast(allPositions, ndim=3)
    dtw[dtw == np.inf] = np.nan
    dtwEuclidean = pd.DataFrame(dtw, columns=names, index=names)

    # If not out file set, write to stdout
    out = sys.stdout if not out else out
    dtwEuclidean.stack().reset_index().to_csv(
        out, header=['id1', 'id2', 'distance'], index=False)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('file',
        help='Simulation info dataset - output of processTUinfo.py')
    custom.add_argument('--out',
        help='File to write Euclidean distances (default: stdout)')
    custom.add_argument('--downsample', type=int, default=None,
        help='File to write Euclidean distances %(default)s')

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
