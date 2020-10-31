#!/usr/bin/env python3

""" Assign ATAC-seq score modifier to each bead. """

import re
import sys
import json
import argparse
import fileinput
from bisect import bisect_right
from utilities import coordinates

def main(infile: str, nbases: int, coordinates: str, scoreColumn: int,) -> None:


    beadDict = buildBead(coordinates['start'], coordinates['end'], nbases)
    with fileinput.input(infile) as fh:
        for record in fh:
            record = record.split()
            chromBed = record[0]
            startBed = int(record[1])
            endBed = int(record[2])
            score = float(record[scoreColumn])
            if inRegion(coordinates['chr'], coordinates['start'],
                        coordinates['end'], chromBed, startBed):
                for pos in range(startBed, endBed):
                    beadStart = findBead(beadDict, pos)
                    if beadStart:
                        beadDict[beadStart].append(score)

    meanBead =  meanBeadScore(beadDict, nbases)
    json.dump(meanBead, sys.stdout)



def splitCoordinates(coordinates):
    """ Split coordinates into ref, start and end """

    ref, start, end = re.split(':|-', a)
    return ref, start, end


def inRegion(chrom, start, end, chromBed, startBed):
    """ Check if BED coordinates are within a genomic region """

    if (chromBed.replace('chr', '') == chrom.replace('chr', '')
            and start <= startBed < end):
        return True
    else:
        return False

def buildBead(start, end, nbases):
    """ Intialise a dictionary of beads with bead start positions  """

    beadDict = {}
    for pos in range(start, end, nbases):
        beadDict[pos] = []
    return beadDict


def findBead(beadDict, pos):
    """ Return the beadStart position of a given base position """

    beadStart = list(beadDict.keys())
    if min(beadStart) > pos > max(beadStart):
        return None
    index =  bisect_right(beadStart, pos) - 1
    return beadStart[index]


def meanBeadScore(beadDict, nbases):
    """ Compute average score for each bead position """

    meanBead = {}
    for i, bead in enumerate(beadDict):
        meanBead[i] = sum(beadDict[bead]) / nbases
    return meanBead


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'coordinates', metavar='CHR:START-END', type=coordinates,
        help='Genomic coordinates to operate on.')
    parser.add_argument(
        'nbases', type=int,
        help='Number of bases to represent 1 bead.')
    parser.add_argument('infile', nargs='?', default=[],
        help='ATAC-seq BED file (default: stdin)')
    parser.add_argument(
        '--scoreColumn', type=int, default=4,
        help='Score column for ATAC seq BED file. (default: %(default)s)')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
