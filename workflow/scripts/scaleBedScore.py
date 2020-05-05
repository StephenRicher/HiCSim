#!/usr/bin/env python3

""" Scale score column of BED file to between 0 and 1 """


import sys
import numpy as np
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(version=__version__, infile=True, in_type='BED')
    parser.set_defaults(function=scaleBED)

    parser.add_argument(
        '--log', default=False, action='store_true',
        help='Log transform scores before scaling')

    return (pct.execute(parser))


def scaleBED(infile, log):

    max_score = maxScore(infile)
    if log:
        max_score = np.log(max_score)
    with pct.open(infile) as f:
        for line in f:
            if line.startswith('#'):
                sys.stdout.write(f'{line}')
                continue
            columns = line.split()
            score = float(columns[4])
            if log:
                score = np.log(score)
            scaled_score = score / max_score
            columns[4] = str(scaled_score)
            line = '\t'.join(columns)
            sys.stdout.write(f'{line}\n')


def maxScore(BED):
    record = 0
    with pct.open(BED) as f:
        for line in f:
            if line.startswith('#'):
                continue
            score = float(line.strip().split()[4])
            if record == 0 or score > max_score:
                max_score = score
            record += 1
    return max_score


if __name__ == '__main__':
    sys.exit(main())
