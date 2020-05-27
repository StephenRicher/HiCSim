#!/usr/bin/env python3

""" Mask a genome with relevant bases """


import sys
import random
import pyCommonTools as pct
from utilities import commaPair
from collections import defaultdict


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__,)
    parser.set_defaults(function=mask)
    parser.add_argument(
        'chromSizes', help='Chromosome sizes file.')
    parser.add_argument(
        '--bed', metavar='BED,CHAR', default=[],
        type=commaPair, action='append',
        help='BED file of regions to mask, paired with masking character.'
        'Call multiple times to add more files.')

    return (pct.execute(parser))


def mask(chromSizes, bed):
    regions = processBed(bed)
    chromosomes = readChromosomes(chromSizes)
    writeMaskedFasta(regions, chromosomes)


def readChromosomes(chrom_sizes):
    chromosomes = {}
    with open(chrom_sizes) as fh:
        for line in fh:
            name, length = line.split()
            chromosomes[name] = int(length)
    return chromosomes


def processBed(beds):
    log = pct.create_logger()
    regions = defaultdict(lambda: defaultdict(list))
    invalid_maskings = ['N', 'B']
    for bed, mask in beds:
        if mask in invalid_maskings:
            log.error(f'Invalid masking character {mask}. '
                      f'Character must not be in {invalid_maskings}.')
            sys.exit(1)
        with open(bed) as fh:
            for line in fh:
                entries = line.split()
                ref = entries[0]
                start = int(entries[1]) + 1 # convert to 1-based
                end = int(entries[2])
                mid = round((start + end)/2)
                regions[ref][mid].append(mask)
    return regions


def writeMaskedFasta(regions, chromosomes):
    for chromosome, length in chromosomes.items():
        sys.stdout.write(f'>{chromosome}\n')
        for index in range(1, length + 1):
            if index in regions[chromosome]:
                base = random.choice(regions[chromosome][index])
            else:
                base = 'N'
            sys.stdout.write(base)
        sys.stdout.write('\n')


if __name__ == '__main__':
    sys.exit(main())
