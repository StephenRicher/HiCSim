#!/usr/bin/env python3

""" Convert custom LAMMPS output to XYZ format """

import sys
import argparse
import fileinput
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def custom2XYZ(file: str):

    with fileinput.input(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if timeStepHeader(line):
                timeStep = next(fh).strip()
            elif nAtomHeader(line):
                nAtoms = next(fh).strip()
            elif coordinateHeader(line):
                print(f'{nAtoms}\nAtoms. Timestep: {timeStep}')
                for line in fh:
                    if timeStepHeader(line):
                        timeStep = next(fh).strip()
                        break
                    else:
                        id, type, x, y, z = line.split()[:5]
                        print(type, x, y, z)


def timeStepHeader(line):
    return line.startswith('ITEM: TIMESTEP')


def coordinateHeader(line):
    return line.startswith('ITEM: ATOMS')


def nAtomHeader(line):
    return line.startswith('ITEM: NUMBER OF ATOMS')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=custom2XYZ)
    parser.add_argument(
        'file', metavar='FILE', nargs='?', default=[],
        help='Lammps custom output (default: stdin)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
