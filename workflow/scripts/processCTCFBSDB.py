#!/usr/bin/env python3

""" Process CTCFBSDB 2.0 prediction tool output to BED """

import os
import re
import sys
import logging
import argparse
import fileinput

__version__ = '1.0.0'

def main(file, threshold=3, **kwargs):
    """ Extract valid CTCFBSDB motifs and convert to BED. """

    with fileinput.input(file) as fh:
        for line in fh:
            motif, seq, name, pos, length, orient, motif_score = line.split()
            if float(motif_score) < threshold:
                continue
            chr, start, end, reps, score = re.sub('[:-]', ' ', name).split()
            start = int(start) + int(pos)
            end = start + int(pos) + int(length)
            sys.stdout.write(f'{chr}\t{start}\t{end}\t{name}\t{score}\t{orient}\n')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='FILE', nargs='?', default=[],
        help='CTCFBSDB prediction output (default: stdin)')
    custom.add_argument(
        '--threshold', default=3, type=int,
        help='Minimum valid motif score (default: %(default)s)')
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
