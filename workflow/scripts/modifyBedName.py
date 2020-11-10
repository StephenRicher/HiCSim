#!/usr/bin/env python3

""" Set name column of BED to first 5 columns to store information """

import sys
import logging
import argparse
import fileinput
from utilities import setDefaults, bedHeader


__version__ = '1.0.0'


def modifyName(file: str):
    """ Read BED file and modify name field """

    with fileinput.input(file) as fh:
        for line in fh:
            if bedHeader(line):
                continue
            chrom, start, end, name, score = line.split()[0:5]
            name = f'{chrom}:{start}-{end}:{name}:{score}'
            print(chrom, start, end, name, score, sep='\t')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='BED file (default: stdin)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(modifyName(**vars(args)))
