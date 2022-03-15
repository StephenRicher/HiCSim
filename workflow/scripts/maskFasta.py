#!/usr/bin/env python3

""" Mask a genome with relevant bases """

import sys
import gzip
import random
import logging
import argparse
import fileinput
from collections import defaultdict
from utilities import setDefaults, createMainParent, coordinates


__version__ = '1.0.0'


def maskFasta(region, beds, nBases: int, seed: float):
    random.seed(seed)
    tracks = processBed(beds, region['chr'])
    bases = []
    n = 0
    for i in range(region['start'], region['end']):
        if i in tracks:
            base = random.choice(tracks[i])
            bases.append(base)
        n += 1
        if n == nBases:
            sys.stdout.write(f'{getBead(bases)}\n')
            bases = []
            n = 0
    if n > 0:
        sys.stdout.write(f'{getBead(bases)}\n')


def processBed(beds, chrom):
    regions = defaultdict(str)
    invalidMaskings = {'N', 'B'}
    for bed in beds:
        with gzip.open(bed, 'rt') as fh:
            for line in fh:
                ref, start, end, mask = line.split()
                if ref != chrom:
                    continue
                mid = round((int(start) + int(end)) / 2)
                regions[mid] += mask
    if set(regions.values()) & invalidMaskings:
        logging.error(f'Invalid masking characters {invalidMaskings}.')
        raise ValueError
    return regions


def getBead(bases):
    if len(bases) == 0:
        bead = 'N'
    else:
        # If both CTCF orientations exist create a bead 'B' to represent both
        if 'F' in bases and 'R' in bases:
            bases = ['B' if base in ['F', 'R'] else base for base in bases]
        bead = random.choice(bases)
    return bead



def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=maskFasta)
    parser.add_argument(
        'region', metavar='CHR:START-END', type=coordinates,
        help='Genomic region (0-based) to operate on.')
    parser.add_argument('beds', nargs='+', help='BED files of regions to mask')
    parser.add_argument(
        '--nBases', default=1000, type=int,
        help='Number of bases to represent 1 bead.')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Seed for random number generator.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
