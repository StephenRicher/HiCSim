#!/usr/bin/env python3

""" Read custom atom group assignments to JSON format """

import sys
import json
import logging
import argparse
import fileinput
from collections import defaultdict

__version__ = '1.0.0'

def main(file : str, **kwargs):

    atomGroups = defaultdict(list)
    with fileinput.input(file) as fh:
        for line in fh:
            if line == 'Atoms\n':
                next(fh) # Skip initial blank line
                atomIdx = 1
                for line in fh:
                    line = line.strip().split()
                    if not line: # Finished Atoms section
                        json.dump(atomGroups, sys.stdout)
                        return 0
                    groups = line[line.index('#') + 1:]
                    for group in groups:
                        atomGroups[group].append(atomIdx)
                    atomIdx += 1


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='DAT', nargs='?', default=[],
        help='Input LAMMPS dat file (default: stdin)')
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
