#!/usr/bin/env python3

""" Simulate a genomic inversion by inverting coordinates and orientation
    between specified position in BED file.
    """

import os
import re
import sys
import logging
import argparse
import fileinput
from utilities import coordinates


__version__ = '1.0.0'


def main(coordinates, file, **kwargs):

    with fileinput.input(file) as fh:
        for line in fh:
            line = line.strip()
            chr, start, end, name, score, orientation = line.split()[:6]
            start = int(start)
            end = int(end)
            if (chr == coordinates['chr']
                    and coordinates['start'] <= start <= coordinates['end']):
                if end > coordinates['end']:
                    logging.warning(
                        f'Entry crosses inversion boundary.\n{line}')
                start2 = coordinates['end'] - (end - coordinates['start'])
                end2 = coordinates['end'] - (start - coordinates['start'])
                start = start2
                end = end2
                if orientation == '+':
                    orientation = '-'
                elif orientation == '-':
                    orientation = '+'
                else:
                    logging.info(f'Entry has no orientation.\n{line}')
            print(chr, start, end, name, score, orientation, sep='\t')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'coordinates', type=coordinates,
        help='Genomic coordinates to invert.')
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
