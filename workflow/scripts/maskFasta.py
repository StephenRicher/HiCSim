#!/usr/bin/env python3

""" Mask a genome with relevant bases """

import sys
import random
import argparse
import logging
import fileinput
from utilities import commaPair, coordinates
from collections import defaultdict


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
                entries = line.split()
                ref = entries[0]
                start = int(entries[1]) + 1 # convert to 1-based
                end = int(entries[2])
                mid = round((start + end)/2)
                regions[ref][mid] += mask

    return regions


def writeMaskedFasta(tracks, chromosomes, region=None):
    for chromosome, length in chromosomes.items():

        if not region:
            indexes = range(1, length + 1)
        elif region['chr'] == chromosome:
            indexes = range(region['start'] + 1, region['end'] + 1)
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


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults()
    custom.add_argument('chromSizes', nargs='?', default=[],
        help='Chromosome sizes file (default: stdin)')
    custom.add_argument(
        '--region', metavar='CHR:START-END', default=None, type=coordinates,
        help='Genomic region (0-based) to operate on.')
    custom.add_argument(
        '--bed', metavar='BED,CHAR', default=[],
        type=commaPair, action='append',
        help='BED file of regions to mask, paired with masking character.'
             'Call multiple times to add more files.')
    custom.add_argument(
        '--seed', default=None, type=float,
        help='Seed for random number generator.')

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
    del args.verbose

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = main(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
