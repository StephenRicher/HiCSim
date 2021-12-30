#!/usr/bin/env python3

""" Read synthetic sequence type and output single bead sequence  """


import re
import sys
import random
import logging
import argparse
import fileinput
import pandas as pd
from argUtils import setDefaults, createMainParent, positiveInt


__version__ = '1.0.0'


def sampleSynthetic(sequence: str, seed: int):
    random.seed(seed)
    validPattern = re.compile(r'^[^:,]+$|^([^:,]+:(\d*\.\d+|\d+),?)+$')
    processedSequence = []
    with fileinput.input(sequence) as fh:
        for i, line in enumerate(fh):
            # Remove all whitespace
            line = re.sub('\s+', '', line)
            if not validPattern.match(line):
                logging.error(f'Invalid patten {line} on line {i}')
                return 1
            probs = parseProbablities(line)
            # Select base weighted by probabilties
            base = random.choices(
                list(probs.keys()), weights=list(probs.values()), k=1)
            processedSequence += base
    print('\n'.join(processedSequence))


def parseProbablities(pattern):
    """ Extract base probabilities from pattern """
    # If no probablities provided set all characters to same probability
    if ':' not in pattern:
        pattern += ':1'
    probs = {}
    for group in pattern.rstrip(',').split(','):
        bases, prob = group.split(':')
        for base in bases:
            if base in probs:
                raise ValueError(f'Duplicate character detected in: {pattern}')
            probs[base] = float(prob)
    return probs


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=sampleSynthetic)
    parser.add_argument(
        'sequence', nargs='?',
        help='Bead sequence to process (default: stdin)')
    parser.add_argument(
        '--seed', default=None, type=positiveInt,
        help='Non-negative integer for seeding random placement '
             'positions (default: %(default)s)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
