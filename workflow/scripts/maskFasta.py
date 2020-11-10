#!/usr/bin/env python3

""" Mask a genome with relevant bases """

import sys
import random
import logging
import argparse
import fileinput
from collections import defaultdict
from utilities import setDefaults, coordinates, bedHeader


__version__ = '1.0.0'


def main(chromSizes, bed, region, seed):

    random.seed(seed)
    tracks = processBed(bed)
    chromosomes = readChromosomes(chromSizes)

    writeMaskedFasta(tracks, chromosomes, region)


def readChromosomes(file):
    chromosomes = {}
    with fileinput.input(file) as fh:
        for line in fh:
            name, length = line.split()
            chromosomes[name] = int(length)
    return chromosomes


def processBed(beds):

    regions = defaultdict(lambda: defaultdict(str))
    invalid_maskings = ['N', 'B']
    for bed, mask in beds:
        if mask in invalid_maskings:
            logging.error(f'Invalid masking character {mask}. '
                          f'Character must not be in {invalid_maskings}.')
            sys.exit(1)
        with open(bed) as fh:
            for line in fh:
                if bedHeader(line): continue
                ref, mid = findMidpoint(line)
                regions[ref][mid] += mask

    return regions


def findMidpoint(line):
    """ Return ref and mid coordinate of BED entry """

    ref, start, end = line.split()[:3]
    mid = round((int(start) + int(end)) / 2)

    return ref, mid


def writeMaskedFasta(tracks, chromosomes, region=None):

    for chromosome, length in chromosomes.items():
        if not region:
            indexes = range(length)
        elif region['chr'] == chromosome:
            indexes = range(region['start'], region['end'])
        else:
            continue

        sys.stdout.write(f'>{chromosome}\n')
        for index in indexes:
            chrTracks = tracks[chromosome]
            if index in chrTracks:
                base = chrTracks[index]
                if len(base) > 1:
                    base = random.choice(base)
            else:
                base = 'N'
            sys.stdout.write(base)
        sys.stdout.write('\n')


def commaPair(value):
    """ Split command seperated pair and return as tuple """

    # Split on last occurence of comma
    bed, maskChar = value.rsplit(',', 1)
    if len(maskChar) != 1:
        raise argparse.ArgumentTypeError(
            f'Masking character {maskChar} must single character.')
    else:
        return (bed, maskChar)


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('chromSizes', nargs='?', default=[],
        help='Chromosome sizes file (default: stdin)')
    parser.add_argument(
        '--region', metavar='CHR:START-END', default=None, type=coordinates,
        help='Genomic region (0-based) to operate on.')
    parser.add_argument(
        '--bed', metavar='BED,CHAR', default=[],
        type=commaPair, action='append',
        help='BED file of regions to mask, paired with masking character.'
             'Call multiple times to add more files.')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Seed for random number generator.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(main(**vars(args)))
