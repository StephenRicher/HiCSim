#!/usr/bin/env python3

""" Calculate Euclidean distance between aligned time series using DTW. """

import sys
import logging
import argparse
import pandas as pd
from tslearn.metrics import cdist_dtw

__version__ = '1.0.0'

def main(file: str, out: str, **kwargs) -> None:

    # Read input and sort by TU then by time and group by TU
    TUinfo = pd.read_csv(file, converters={'pos': eval})
    TUnames = TUinfo['TU'].unique()
    TUinfo = TUinfo.sort_values(['TU', 'time']).groupby(['TU'])

    # Loop through TUs and add to timeseries data structure
    allPositions = []
    for name, TUtimeseries in TUinfo['pos']:
        allPositions.append(TUtimeseries.to_list())

    # Compute time aligned Euclidean distance

    dtwEuclidean = pd.DataFrame(
        cdist_dtw(allPositions), columns=TUnames, index=TUnames)

    # If not out file set, write to stdout
    out = sys.stdout if not out else out
    dtwEuclidean.stack().reset_index().to_csv(
        out, header=['TU1', 'TU2', 'distance'], index=False)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('file', help='TU dataset - output of processTUinfo.py')
    custom.add_argument('--out',
        help='File to write Euclidean distances (default: stdout)')
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
