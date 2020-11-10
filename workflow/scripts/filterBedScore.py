#!/usr/bin/env python3

""" Sample BED file based on scaled SCORE """

import sys
import random
import argparse
import fileinput
from utilities import setDefaults, bedHeader


__version__ = '1.0.0'


def filterBed(file: str, seed: float):

    random.seed(seed)
    with fileinput.input(file) as fh:
        for line in fh:
            if bedHeader(line): continue
            score = float(line.split()[4])
            assert 0 <= score <= 1, "Score not between 0 and 1"
            if random.random() < score:
                print(line)


def parseArgs():
    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='Input BED file (default: stdin)')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Initialize the random number generator (default: %(default)s)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(filterBed(**vars(args)))
