#!/usr/bin/env python3

""" Split BED file records by record orientation. """

import sys
import logging
import argparse
import fileinput
from contextlib import ExitStack
from utilities import setDefaults, bedHeader


__version__ = '1.0.0'


def main(file: str, forward: str, reverse: str, minRep: int):

    with ExitStack() as stack:

        fh = stack.enter_context(fileinput.input(file))
        outForward = stack.enter_context(open(forward, 'w'))
        outReverse = stack.enter_context(open(reverse, 'w'))

        for i, line in enumerate(fh):
            if bedHeader(line):
                continue
            columns = line.split()
            rep = int(columns[3])
            if rep < minRep:
                continue
            orientation = columns[5]
            if orientation == '+':
                out = outForward
            elif oritentation == '-':
                out = outReverse
            else:
                logging.info(f'No orientation on line {i+1} - skipping.')
                continue

            print(line, end='', file=out)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    custom.add_argument(
        'file', metavar='BED', nargs='?', default=[],
        help='Input BED file (default: stdin)')
    custom.add_argument(
        '--minRep', type=int, default=1,
        help='Minimum replicates required to write (default: %(default)s)')
    requiredNamed = custom.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--reverse', required=True,
        help='Reverse orientation BED output.')
    requiredNamed.add_argument(
        '--forward', required=True,
        help='Forward orientation BED output.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(modifyName(**vars(args)))
