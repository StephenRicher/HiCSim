#!/usr/bin/env python3

""" Process CTCFBSDB 2.0 prediction tool output to BED """

import os
import re
import sys
import logging
import argparse
import fileinput
from utilities import setDefaults, createMainParent

__version__ = '1.0.0'


def main(file, threshold=3, **kwargs):
    """ Extract valid CTCFBSDB motifs and convert to BED. """

    with fileinput.input(file) as fh:
        for line in fh:
            motif, seq, name, pos, length, orient, motif_score = line.split()
            if float(motif_score) < threshold:
                continue
            chr, start, end, reps, score = re.sub('[:-]', ' ', name).split()
            start = int(start) + int(pos)
            end = start + int(pos) + int(length)
            sys.stdout.write(f'{chr}\t{start}\t{end}\t{name}\t{score}\t{orient}\n')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=main)
    parser.add_argument(
        'file', metavar='FILE', nargs='?', default=[],
        help='CTCFBSDB prediction output (default: stdin)')
    parser.add_argument(
        '--threshold', default=3, type=int,
        help='Minimum valid motif score (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
