#!/usr/bin/env python3


""" Sample BED file based on scaled SCORE """

import sys
import random
import logging
import argparse
import fileinput

__version__ = '1.0.0'


def main(file : str, seed : float, **kwarfs):

    random.seed(seed)

    with fileinput.input(file) as fh:
        for line in fh:
            if line.startswith('#'):
                sys.stdout.write(f'{line}')
                continue
            columns = line.split()
            score = float(columns[4])
            if random.random() < score:
                columns[4] = str(score)
                line = '\t'.join(columns)
                sys.stdout.write(f'{line}\n')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='Input BED file (default: stdin)')
    custom.add_argument(
        '--seed', default=None, type=float,
        help='Initialize the random number generator (default: %(default)s)')
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
