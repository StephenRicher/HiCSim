#!/usr/bin/env python3

""" Read custom atom group assignments to JSON format """

import sys
import json
import argparse
import fileinput
from utilities import setDefaults
from collections import defaultdict


__version__ = '1.0.0'


def getAtomGroups(file: str):

    atomGroups = defaultdict(list)
    with fileinput.input(file) as fh:
        for line in fh:
            if line == 'Atoms\n':
                next(fh) # Skip initial blank line
                atomIdx = 1
                for line in fh:
                    line = line.strip().split()
                    if not line: # Finished Atoms section
                        json.dump(atomGroups, sys.stdout)
                        return 0
                    groups = line[line.index('#') + 1:]
                    if groups:
                        for group in groups:
                            # Group all CTCF sites by direction
                            if group.startswith(('F-', 'R-', 'B-')):
                                atomGroups[group[0]].append(atomIdx)
                            atomGroups[group].append(atomIdx)
                    else:
                        atomGroups["None"].append(atomIdx)
                    atomIdx += 1


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'file', metavar='DAT', nargs='?', default=[],
        help='Input LAMMPS dat file (default: stdin)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(getAtomGroups(**vars(args)))
