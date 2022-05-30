#!/usr/bin/env python3

""" Sample BED file based on scaled SCORE """

import sys
import random
import logging
import argparse
import fileinput
from utilities import setDefaults, createMainParent, bedHeader


__version__ = '1.0.0'


def filterBed(bed: str, char: list, seed: float, retentionProb: float):

    random.seed(seed)
    assert len(char) <= 2, "Maximum 2 characters."
    assert 0 <= retentionProb <= 1, 'Filter probability not between 0 and 1.'
    pairChar = len(char) == 2
    with open(bed) as fh:
        for line in fh:
            if bedHeader(line):
                continue
            columns = line.strip().split()
            if random.random() < retentionProb:
                if (pairChar) and (len(columns) >=6) and (columns[5] == '-'):
                    c = char[1]
                else:
                    c = char[0]
                print(columns[0], columns[1], columns[2], c)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=filterBed)
    parser.add_argument(
        'bed', metavar='BED', help='Input BED file')
    parser.add_argument(
        'char', nargs='+',
        help='Masking character for BED intervals. If 2 provided, they will '
             'be used for forward / reverse orientation respectively.')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Initialize the random number generator (default: %(default)s)')
    parser.add_argument(
        '--retentionProb', default=0.5, type=float,
        help='Set retention probability for each'
             'interval (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
