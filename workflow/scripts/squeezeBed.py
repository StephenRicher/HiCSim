#!/usr/bin/env python3

""" Squeeze/expand intervals of BED file to specified length """

import sys
import math
import random
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(version=__version__)
    parser.set_defaults(function=squeezeBED)

    parser.add_argument(
        'length', type=int,
        help='Initialize the random number generator.')
    parser.add_argument(
        'bed', metavar='BED', nargs='?',
        help='Input BED file.')

    return (pct.execute(parser))


def squeezeBED(bed, length):

    with pct.open(bed) as f:
        for line in f:
            if line.startswith('#'):
                sys.stdout.write(f'{line}')
                continue
            entry = line.strip().split()
            interval_midpoint = round((int(entry[2]) + int(entry[1])) / 2)
            above_mid = round(int(length / 2))
            below_mid = length - above_mid
            entry[1] = str(interval_midpoint - below_mid)
            entry[2] = str(interval_midpoint + above_mid)
            print('\t'.join(entry))


if __name__ == '__main__':
    sys.exit(main())
