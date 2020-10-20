#!/usr/bin/env python3

""" Compute mean and standard deviation of DTW distance for a given seperation  """

import sys
import logging
import argparse
import pandas as pd
from typing import List


__version__ = '1.0.0'


def main(files: List, out: str, **kwargs) -> None:

    normalise = pd.concat((pd.read_csv(file) for file in files))
    normalise['seperation'] = abs(normalise['id1'] - normalise['id2'])
    normalise = (normalise
        .set_index(['id1', 'id2'])
        .groupby('seperation').agg(['mean', 'std'])
        .reset_index())
    normalise.columns = ['seperation', 'mean', 'std']
    normalise.to_csv(out, index=False)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'files', nargs='+', help='Input DTW distances.')
    custom.add_argument(
        '--out', default=sys.stdout,
        help='Output file for normalisation factors (default: stdout')
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
