#!/usr/bin/env python3

""" Subset ATAC JSON bead modifier by genomic region """

import sys
import json
import logging
import argparse
from utilities import getBead, coordinates, readJSON


def main(infile: str, coordinates: dict, nbases: int) -> None:

    if not validCoordinates(coordinates['start'], coordinates['end'], nbases):
        logging.error('Error not valid coordinates.')
        return 1

    beadDict = readJSON(infile)

    chrom = re.sub('^chr', '', coordinates['chr'])
    startBead = getBead(coordinates['start'], nbases)
    endBead = getBead(coordinates['end'], nbases)
    subsetBead = {}

    for bead in range(startBead, endBead):
        try:
            subsetBead[bead] = beadDict[chrom][bead]
        except KeyError:
            subsetBead[bead] = 0

    json.dump(subsetBead, sys.stdout)


def validCoordinates(start, end, nbases):
    """ Check if start/end positions are multiples of nbases """
    return start % nbases == end % nbases == 0


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        '--coordinates', required=True, metavar='CHR:START-END', type=coordinates,
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
