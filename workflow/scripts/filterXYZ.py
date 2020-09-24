#!/usr/bin/env python3

""" Extract specific atom group from XYZ file """

import os
import re
import sys
import json
import logging
import argparse
import fileinput

__version__ = '1.0.0'


def main(group, atomGroups, file, **kwargs):

    atomGroupsDict = readJSON(atomGroups)
    try:
        groupIdxs = atomGroupsDict[group]
    except KeyError:
        logging.error(f'{group} group not present in {atomGroups}.')
        sys.exit(1)

    nAtoms = len(groupIdxs)
    with fileinput.input(file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if coordinateLine(line):
                if counter in groupIdxs:
                    print(line)
                counter += 1
            elif nAtomLine(line):
                print(nAtoms)
            elif commentLine(line):
                counter = 1
                print(line)


def commentLine(line):
    return line.startswith('Atoms')


def coordinateLine(line):
    return len(line.split()) == 4


def nAtomLine(line):
    try:
        atoms = int(line.strip())
        return True
    except ValueError:
        return False


def readJSON(file):
    """ Read JSON encoded data to dictionary """
    with open(file) as fh:
        return json.load(fh)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'group',
        help='Desired group to filter.')
    custom.add_argument(
        'atomGroups',
        help='Atom group assignments in JSON format.')
    custom.add_argument(
        'file', metavar='XYZ', nargs='?', default=[],
        help='Coordinate file generated by lammps (default: stdin)')
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
