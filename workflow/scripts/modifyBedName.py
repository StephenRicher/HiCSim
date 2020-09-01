#!/usr/bin/env python3

""" Set name column of BED to first 5 columns to store information """

import sys
import logging
import argparse
import fileinput

__version__ = '1.0.0'


def main(file  : str, **kwargs):
    """ Read BED file and modify name field """

    with fileinput.input(file) as fh:
        for line in fh:
            chrom, start, end, name, score = line.strip('\n').split('\t')[0:5]
            name = f'{chrom}:{start}-{end}:{name}:{score}'
            print(chrom, start, end, name, score, sep='\t')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='BED file (default: stdin)')
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
