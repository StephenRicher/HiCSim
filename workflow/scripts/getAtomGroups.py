#!/usr/bin/env python3

""" Read beadID to type assignments of polymer """

import sys
import json
import argparse
import fileinput
from typing import List
from collections import defaultdict
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def getAtomGroups(beads: str, nMonomers: int, TUs: List):

    atomGroups = defaultdict(list)
    beadID = 1
    with fileinput.input(beads) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # Every bead is a DNA bead
            atomGroups['DNA'].append(beadID)
            atomGroups[line].append(beadID)
            if line in TUs:
                atomGroups['TU'].append(beadID)
            beadID += 1
    for monomer in range(nMonomers):
        atomGroups['TF'].append(beadID)
        beadID += 1

    json.dump(atomGroups, sys.stdout)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=getAtomGroups)
    parser.add_argument(
        'beads', nargs='?', default=[], help='Beads file. (default: stdin)')
    parser.add_argument(
        '--nMonomers', default=0, type=int,
        help='Set number of monomers (TFs) for simulation.')
    parser.add_argument(
        '--TUs', nargs='*', default=[],
        help='Bead types that should be treated as transcriptional units.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
