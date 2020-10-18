#!/usr/bin/env python3

""" Calculate TU-TU correlation """

import sys
import logging
import argparse
import pandas as pd
from scipy.stats import pearsonr


__version__ = '1.0.0'

def main(file: str, **kwargs) -> None:

    allTUinfo = pd.read_csv(file)
    allTUinfo = allTUinfo.pivot(index='time', columns='id', values='active')

    # Compute Pearson correlation coefficient and convert to long format
    correlations = allTUinfo.corr(method='pearson').stack()
    correlations.index.names = ['row', 'column']
    correlations = correlations.reset_index()
    correlations.columns = ['row', 'column', 'r']

    # Combpute corresponding p-values and convert to seperate dataframe
    pvalues = allTUinfo.corr(method=pearsonr_pval).stack()
    pvalues.index.names = ['row', 'column']
    pvalues = pvalues.reset_index()
    pvalues.columns = ['row', 'column', 'p-value']

    # Combine to include p-values
    pd.merge(correlations, pvalues).to_csv(sys.stdout, index=False)


def pearsonr_pval(x,y):
    return pearsonr(x,y)[1]


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('file', help='TU dataset - output of processTUinfo.py')
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
