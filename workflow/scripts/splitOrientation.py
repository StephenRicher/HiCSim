#!/usr/bin/env python3

""" Split BED file records by record orientation. """

import sys
import logging
import argparse
import fileinput
from contextlib import ExitStack
from utilities import setDefaults, bedHeader

__version__ = '1.0.0'


def splitStrand(file: str, forward: str, reverse: str):

    with ExitStack() as stack:
        fh = stack.enter_context(fileinput.input(file))
        outForward = stack.enter_context(open(forward, 'w'))
        outReverse = stack.enter_context(open(reverse, 'w'))
        for i, line in enumerate(fh):
            if bedHeader(line):
                continue
            chrom, start, end, name, score, strand = line.split()[:6]
            if strand == '+':
                out = outForward
            elif strand == '-':
                out = outReverse
            else:
                logging.info(f'No orientation on line {i+1} - skipping.')
                continue
            print(chrom, start, end, name, score, strand, sep='\t', file=out)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='Input BED file (default: stdin)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--reverse', required=True,
        help='Reverse strand BED output.')
    requiredNamed.add_argument(
        '--forward', required=True,
        help='Forward strand BED output.')

    return setDefaults(parser)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(splitStrand(**vars(args)))
