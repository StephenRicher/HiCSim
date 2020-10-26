#!/usr/bin/env python3

""" Merge replicate Pandas dataframes """

import sys
import argparse
import pandas as pd
from typing import List


def main(infiles: List, out: str) -> None:

    mergedDF = []
    for rep, file in enumerate(infiles):
        df = pd.read_csv(file)
        df['rep'] = rep
        mergedDF.append(df)

    pd.concat(mergedDF).to_csv(out, index=False)


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infiles', nargs='*',
        help='Input CSV files to merge.')
    parser.add_argument('--out', default=sys.stdout,
        help='Merged output file.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
