#!/usr/bin/env python3


""" Sample BED file based on scaled SCORE """

import sys
import random
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(version=__version__, infile=True, in_type='BED')
    parser.set_defaults(function=filterBED)


    parser.add_argument(
        '--seed', default=None, type=int,
        help='Initialize the random number generator.')

    return (pct.execute(parser))


def filterBED(infile, seed):

    random.seed(seed)

    with pct.open(infile) as f:
        for line in f:
            if line.startswith('#'):
                sys.stdout.write(f'{line}')
                continue
            columns = line.split()
            score = float(columns[4])
            if random.random() < score:
                columns[4] = str(score)
                line = '\t'.join(columns)
                sys.stdout.write(f'{line}\n')


if __name__ == '__main__':
    sys.exit(main())
