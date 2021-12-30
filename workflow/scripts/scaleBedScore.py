#!/usr/bin/env python3

""" Scale score column of BED file to between 0 and 1 """

import sys
import logging
import argparse
import numpy as np
from utilities import setDefaults, transformScore, bedHeader

__version__ = '1.0.0'


def scaleBed(file: str, transform: str):

    maxScore = getMaxScore(file)
    maxScore = transformScore(maxScore, transform)

    with open(file) as fh:
        for line in fh:
            if bedHeader(line): continue
            columns = line.strip().split()
            score = float(columns[4])
            score = transformScore(score, transform)
            score = score / maxScore
            columns[4] = str(score)
            line = '\t'.join(columns)
            print(line)


def getMaxScore(file):
    """ Determine maxScore """

    with open(file) as fh:
        for line in fh:
            if bedHeader(line): continue
            try:
                score = float(line.split()[4])
                maxScore = max(score, maxScore)
            except UnboundLocalError:
                maxScore = score
    return maxScore


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('file', help='Input BED file')
    parser.add_argument(
        '--transform', default='none', choices=['log', 'sqrt', 'none'],
        help='Transform to apply to BED scores (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(scaleBed(**vars(args)))
