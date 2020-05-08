#!/usr/bin/env python3

""" Remove unbonded monomers from end of XYZ file before making contact map """

import sys
import argparse
import pyCommonTools as pct
from timeit import default_timer as timer


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'xyz', metavar='XYZ',
        help='Coordinate file generated by lammps.')
    parser.add_argument(
        'nMonomers', type=int,
        help='Number of monomer beads to filter from end of XYZ file')

    parser.set_defaults(function=filterXYZ)

    return (pct.execute(parser))


def filterXYZ(xyz, nMonomers):
    nMonomers = 200
    with open(xyz) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if nAtomLine(line):
                nBound = int(line) - nMonomers
                print(nBound)
            elif commentLine(line):
                counter = 0
                print(line)
            elif coordinateLine(line) and counter < nBound:
                print(line)
                counter += 1


def commentLine(line):
    return line.startswith('Atoms')


def coordinateLine(line):
    return len(line.split()) == 4


def nAtomLine(line):
    try:
        int(line.strip())
        return True
    except ValueError:
        return False


if __name__ == '__main__':
    log = pct.create_logger()
    start = timer()
    RC = main()
    end = timer()
    log.info(f'Total time elapsed: {end - start} seconds.')
    sys.exit(RC)
