#!/usr/bin/env python3

""" Script to read a single record FASTA file and compress N bases to 1 bead """

import sys
import random
import logging
import argparse
import fileinput

__version__ = '1.0.0'


def main(file : str, nbases : int, **kwargs):

    with fileinput.input(file) as fh:
        bases = []
        for i, line in enumerate(fh):
            if line.startswith('>'):
                if i != 0:
                    logging.warning('Second FASTA record detected - skipping.')
                    break
            else:
                if i == 0:
                    logging.error('FASTA does not begin with ">".')
                    sys.exit(1)
                else:
                    for base in line.strip():
                        bases.append(base)
                        if len(bases) == nbases:
                            sys.stdout.write(f'{get_bead(bases)}\n')
                            bases = []
        if len(bases) > 0:
            sys.stdout.write(f'{get_bead(bases)}\n')


def get_bead(bases):
    # Remove any 'N' bases from list
    bases = [base for base in bases if base != 'N']
    if len(bases) == 0:
        bead = 'N'
    else:
        # If both CTCF orientations exist create a bead 'B' to represent both
        if 'F' in bases and 'R' in bases:
            bases = ['B' if base in ['F', 'R'] else base for base in bases]
        bead = random.choice(bases)
    return bead


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='FASTA', nargs='?', default=[],
        help='Input FASTA file (default: stdin)')
    custom.add_argument(
        '--nbases', default=1000, type=int,
        help='Number of bases to represent 1 bead.')
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
