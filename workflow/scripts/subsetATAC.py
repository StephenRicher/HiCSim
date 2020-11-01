#!/usr/bin/env python3

""" Subset ATAC JSON bead modifier by genomic region """

import re
import sys
import json
import logging
import argparse
from utilities import getBead, coordinates, readJSON


def main(infile: str, region: dict, nbases: int) -> None:

    if not validCoordinates(region['start'], region['end'], nbases):
        logging.error('Error not valid coordinates.')
        return 1

    beadDict = readJSON(infile)

    chrom = re.sub('^chr', '', region['chr'])
    startBead = getBead(region['start'], nbases)
    endBead = getBead(region['end'], nbases)
    subsetBead = {}

    for i, bead in enumerate(range(startBead, endBead)):
        try:
            subsetBead[i] = beadDict[chrom][bead]
        except KeyError:
            subsetBead[i] = 0

    json.dump(subsetBead, sys.stdout)


def validCoordinates(start, end, nbases):
    """ Check if start/end positions are multiples of nbases """
    return start % nbases == end % nbases == 0


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        '--region', required=True, metavar='CHR:START-END', type=coordinates,
        help='Genomic coordinates to operate on.')
    parser.add_argument(
        '--nbases', required=True, type=int,
        help='Number of bases to represent 1 bead.')
    parser.add_argument('infile', nargs='?', default=[],
        help='ATAC JSON bead file (default: stdin)')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
