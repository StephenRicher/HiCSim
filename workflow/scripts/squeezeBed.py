#!/usr/bin/env python3

""" Reform intervals of BED file to specified length """

import sys
import math
import logging
import argparse
import fileinput

__version__ = '1.0.0'


def main(file : str, length : int, **kwargs):

    with fileinput.input(file) as fh:
        for line in fh:
            if line.startswith('#'):
                sys.stdout.write(f'{line}')
                continue
            entry = line.strip().split()
            interval_midpoint = round((int(entry[2]) + int(entry[1])) / 2)
            above_mid = round(int(length / 2))
            below_mid = length - above_mid
            entry[1] = str(interval_midpoint - below_mid)
            entry[2] = str(interval_midpoint + above_mid)
            print('\t'.join(entry))


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'length', type=int,
        help='BED interval size to reform to.')
    custom.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='Input BED file (default: stdin)')
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
