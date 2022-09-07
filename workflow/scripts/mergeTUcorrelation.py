#!/usr/bin/env python3

""" Merge TU correlation data """


import sys
import argparse
import pandas as pd
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import fdrcorrection
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def mergeTUcorrelation(infiles: list, out: str):
    df = []
    for file in infiles:
        df.appennd(pd.read_csv(file))
    df = pd.concat(df)

    merged = df.groupby(['row', 'column'])['p'].apply(
        lambda x: combine_pvalues(x)[1]).reset_index()
    merged['p(adj)'] = fdrcorrection(merged['p'])[1]
    merged['count'] = df.groupby(['row', 'column']).size().values
    merged['r'] = df.groupby(['row', 'column'])['r'].median().values
    merged.to_csv(out, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=mergeTUcorrelation)
    parser.add_argument('infiles', nargs='*', help='Input files to merge')
    parser.add_argument(
        '--out', default=sys.stdout, help='Merged output (default: stdout)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
