#!/usr/bin/env python3

""" Convert custom LAMMPS dump output to XYZ format """

import os
import re
import sys
import logging
import argparse
import fileinput
import collections

__version__ = '1.0.0'


def main(file, **kwargs):

    with fileinput.input(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if timeStepHeader(line):
                timeStep = int(next(fh).strip())
            elif nAtomHeader(line):
                nAtoms = int(next(fh).strip())
            elif coordinateHeader(line):
                print(f'{nAtoms}\nAtoms. Timestep: {timeStep}')
                for line in fh:
                    try:
                        id, type, x, y, z, ix, iy, iz = line.strip().split()
                        print(type, x, y, z)
                    # Continue reading coordinates until next timeStepHeader
                    except ValueError:
                        timeStep = int(next(fh).strip())
                        break



def timeStepHeader(line):
    return line.startswith('ITEM: TIMESTEP')


def coordinateHeader(line):
    return line.startswith('ITEM: ATOMS')


def nAtomHeader(line):
    return line.startswith('ITEM: NUMBER OF ATOMS')


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='FILE', nargs='?', default=[],
        help='Lammps custom dump outputRemoteTowe (default: stdin)')
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
