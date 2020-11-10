#!/usr/bin/env python3

""" Scale score column of BED file to between 0 and 1 """

import sys
import logging
import argparse
import numpy as np
from utilities import setDefaults, transformScore, bedHeader

__version__ = '1.0.0'


def scaleBed(file : str, transform : str, scoreColumn: int):

    maxScore = getMaxScore(file, scoreColumn)
    maxScore = transformScore(maxScore, transform)

    with open(file) as fh:
        for line in fh:
            if bedHeader(line): continue
            columns = line.split()
            score = float(columns[scoreColumn])
            score = transformScore(score, transform)
            scaledScore = score / maxScore
            columns[scoreColumn] = str(scaledScore)
            line = '\t'.join(columns)
            sys.stdout.write(f'{line}\n')


def getMaxScore(file, scoreColumn=4):
    """ Determine maxScore """

    with open(file) as fh:
        for line in fh:
            if bedHeader(line): continue
            try:
                score = float(line.split()[scoreColumn])
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
    parser.add_argument(
        '--scoreColumn', type=int, default=4,
        help='Specify score column to scale (default: %(default)s)')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(scaleBed(**vars(args)))
