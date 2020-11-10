#!/usr/bin/env python3

""" Subset ATAC JSON bead modifier by genomic region """

import re
import sys
import json
import logging
import argparse
from utilities import setDefaults, getBead, coordinates, readJSON

__version__ = '1.0.0'


def subsetATAC(infile: str, region: dict, nBases: int) -> None:

    if not validCoordinates(region['start'], region['end'], nBases):
        logging.error('Error not valid coordinates.')
        return 1

    beadDict = readJSON(infile)

    chrom = re.sub('^chr', '', region['chr'])
    startBead = getBead(region['start'], nBases)
    endBead = getBead(region['end'], nBases)
    
    subsetBead = {}
    for i, bead in enumerate(range(startBead, endBead)):
        try:
            subsetBead[i] = beadDict[chrom][str(bead)]
        except KeyError:
            subsetBead[i] = 0

    json.dump(subsetBead, sys.stdout)


def validCoordinates(start: int, end: int, nBases: int) -> bool:
    """ Check if start/end positions are multiples of nBases """
    return start % nBases == end % nBases == 0


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infile', nargs='?', default=[],
        help='ATAC JSON bead file (default: stdin)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--region', metavar='CHR:START-END', required=True, type=coordinates,
        help='Genomic coordinates to operate on.')
    requiredNamed.add_argument(
        '--nBases', required=True, type=int,
        help='Number of bases to represent 1 bead.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(subsetATAC(**vars(args)))
