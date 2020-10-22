#!/usr/bin/env python3

""" Normalise DTW distances for specific seperation with normalisation factors  """

import sys
import logging
import argparse
import pandas as pd
from typing import List
import scipy.stats as st

__version__ = '1.0.0'


def main(normalisationFactors: str, dtwDistances: str, out: str, **kwargs) -> None:

    normaliseFactors = pd.read_csv(normalisationFactors)
    dtw = pd.read_csv(dtwDistances)

    dtw['seperation'] = abs(dtw['id1'] - dtw['id2'])
    dtw = pd.merge(dtw, normaliseFactors, on='seperation')
    # Multiple by -1 so that positive Z score are for closer distances
    dtw['z'] = ((dtw['distance'] - dtw['mean']) / dtw['std']) * -1

    # Rescale distance to contact probability
    dtw['distance'] = st.norm.cdf(dtw['z'])

    dtw.to_csv(out, index=False)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('normalisationFactors',
        help='Input noramlisation factors.')
    custom.add_argument('dtwDistances',
        help='DTW distances to normalise.')
    custom.add_argument(
        '--out', default=sys.stdout,
        help='Output file for normalisation factors (default: stdout)')
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
