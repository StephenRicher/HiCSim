#!/usr/bin/env python3

""" Average contact frequency matrices and plot heatmap """

import sys
import json
import logging
import pandas as pd
import argparse
from typing import List
from utilities import readJSON

__version__ = '1.0.0'


def main(files: List, **kwargs) -> None:

    allBeadDistribution = {}
    for i, file in enumerate(files):
        beads = readJSON(file)
        # Set all beads in first iteration.
        if i == 0:
            for atom in beads['DNA']:
                allBeadDistribution[atom] = 0
        for TU in beads['TU']:
            allBeadDistribution[TU] += 1

    allBeadDistribution = pd.DataFrame(allBeadDistribution, index=[0])
    allBeadDistribution.to_csv(sys.stdout, index=False)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('files', nargs='+', help='Atom JSON file.')
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
