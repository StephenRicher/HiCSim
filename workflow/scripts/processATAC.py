#!/usr/bin/env python3

""" Assign ATAC-seq score modifier to each bead. """

import re
import sys
import json
import argparse
import fileinput
import numpy as np
from collections import defaultdict
from utilities import setDefaults, getBead, transformScore


def processATAC(
        infile: str, transform: str, nBases: int,
        percentile: float, precision: int):

    beadDict = defaultdict(lambda: defaultdict(float))
    with fileinput.input(infile) as fh:
        for record in fh:
            chrom, start, end, score = splitBedgraph(record, strip=True)
            if score == 0:
                continue
            # If start/end map to same bead - add full score
            if getBead(start, nbases) == getBead(end, nbases):
                beadDict[chrom][getBead(start, nbases)] += score
            # Else loop through each base and add score individually
            else:
                length = end - start
                for pos in range(start, end):
                    bead = getBead(pos, nbases)
                    beadDict[chrom][bead] += (score / length)

    maxBeadScore = percentileScore(beadDict, percentile)
    maxBeadScore = transformScore(maxBeadScore, transform)

    for chrom, beads in beadDict.items():
        for bead, score in beads.items():
            score = transformScore(score, transform)
            # Ensure score is not greater than 1
            scaledScore = min(1, score / maxBeadScore)
            # Roud to n decimal places
            if precision:
                scaledScore = round(scaledScor, precision)
            beadDict[chrom][bead] = scaledScore

    json.dump(beadDict, sys.stdout)


def splitBedgraph(record, strip=False):
    """ Return bedgraph values correctly typed """

    chrom, start, end, score = record.split()[0:4]
    if strip:
        chrom = re.sub('^chr', '', chrom)
    start = int(start)
    end = int(start)
    score = float(score)
    return chrom, start, end, score


def percentileScore(beadDict, q=99):
    """ Return percentile score of all score in beadDict """

    allScores = [list(d.values()) for d in beadDict.values()]
    return np.percentile(np.concatenate(allScores), q)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infile', nargs='?', default=[],
        help='ATAC-seq bedgraph file (default: stdin)')
    parser.add_argument(
        '--transform', default='none', choices=['none', 'log', 'sqrt'],
        help='Transform to apply to scores (default: %(default)s)')
    parser.add_argument(
        '--percentile', type=float, default=99,
        help='Percentile score to scale by (default: %(default)s)')
    parser.add_argument(
        '--precision', type=int,
        help='Precision to round modifier score (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--nBases', required=True, type=int,
        help='Number of bases to represent 1 bead.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(processATAC(**vars(args)))
