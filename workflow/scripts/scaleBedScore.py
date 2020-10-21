#!/usr/bin/env python3

""" Scale score column of BED file to between 0 and 1 """

import sys
import logging
import argparse
import fileinput
import numpy as np

__version__ = '1.0.0'


def main(file : str, transform : str, **kwargs):

    max_score = maxScore(file)
    if transform == 'sqrt':
        max_score = np.sqrt(max_score)
    elif transform == 'log':
        max_score = np.log(max_score)

    with fileinput.input(file) as fh:
        for line in fh:
            if line.startswith('#'):
                sys.stdout.write(f'{line}')
                continue
            score = float(line.split()[4])
            if transform == 'sqrt':
                score = np.sqrt(score)
            elif transform == 'log':
                score = np.log(score)
            scaled_score = score / max_score
            columns[4] = str(scaled_score)
            line = '\t'.join(columns)
            sys.stdout.write(f'{line}\n')


def maxScore(BED):
    record = 0
    with fileinput.input(BED) as f:
        for line in f:
            if line.startswith('#'):
                continue
            score = float(line.strip().split()[4])
            if record == 0 or score > max_score:
                max_score = score
            record += 1
    return max_score


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='Input BED file (default: stdin)')
    custom.add_argument(
        '--transform', default='none', choices=['log', 'sqrt', 'none'],
        help='Transform to apply to BED scores (default: %(default)s)')
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
