#!/usr/bin/env python3

""" Split BED file records by record orientation. """

import sys
import logging
import argparse
import fileinput
from contextlib import ExitStack

__version__ = '1.0.0'

def main(file : str, forward :str, reverse : str, min_rep : int, **kwargs):

    with ExitStack() as stack:

        fh = stack.enter_context(fileinput.input(file))
        out_forward = stack.enter_context(open(forward, 'w'))
        out_reverse = stack.enter_context(open(reverse, 'w'))

        for line in fh:
            columns = line.strip().split()
            rep = int(columns[3])
            if rep < min_rep:
                continue
            orientation = columns[5]
            if orientation == '+':
                out = out_forward
            else:
                out = out_reverse

            print(line, end='', file = out)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='Input BED file (default: stdin)')
    custom.add_argument(
        '--min_rep', type=int, default=1,
        help='Minimum number of replicates required per '
             'record to write (default: %(default)s)')
    requiredNamed = custom.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--reverse', required=True,
        help='Output file for reverse orientation intervals.')
    requiredNamed.add_argument(
        '--forward', required=True,
        help='Output file for forward orientation intervals.')
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
