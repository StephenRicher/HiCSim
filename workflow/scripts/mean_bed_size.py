#!/usr/bin/env python3

""" Calculate median interval size of BED file. """

import sys
import argparse
import pyCommonTools as pct
from statistics import median
from timeit import default_timer as timer


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True, in_type='BED',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.set_defaults(function=median_bed)

    return (pct.execute(parser))


def median_bed(infile):

    log = pct.create_logger()

    interval_lengths = []
    with pct.open(infile) as f:
        for line in f:
            line = line.strip().split()
            start = int(line[1])
            end = int(line[2])
            length = end - start

            interval_lengths.append(length)
    print(interval_lengths)

    print(median(interval_lengths))


if __name__ == '__main__':
    log = pct.create_logger()
    start = timer()
    RC = main()
    end = timer()
    log.info(f'Total time elapsed: {end - start} seconds.')
    sys.exit(RC)
