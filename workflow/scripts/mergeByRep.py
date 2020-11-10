#!/usr/bin/env python3

""" Merge pandas dataframes """


import sys
import argparse
import pandas as pd
from typing import List
from utilities import setDefaults


__version__ = '1.0.0'


def mergeByRep(infiles: List, out: str):

    mergedDF = []
    for rep, file in enumerate(infiles):
        df = pd.read_csv(file)
        df['rep'] = rep
        mergedDF.append(df)

    pd.concat(mergedDF).to_csv(out, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infiles', nargs='*', help='Input files to merge')
    parser.add_argument(
        '--out', default=sys.stdout, help='Merged output (default: stdout)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(mergeByRep(**vars(args)))
