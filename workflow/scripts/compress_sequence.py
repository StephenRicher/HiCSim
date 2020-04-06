#!/usr/bin/env python3

""" Script to read a single record FASTA file and compress N bases to 1 bead """

import sys
import argparse
import pyCommonTools as pct
from collections import Counter
from timeit import default_timer as timer


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True, in_type='FASTA',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-n', '--nbases', default=100, type=int,
        help='Number of base to represent 1 bead.')
    parser.set_defaults(function=compress)

    return (pct.execute(parser))


def compress(infile, nbases):

    log = pct.create_logger()

    with pct.open(infile) as f:

        bases = Counter()
        i = 0
        for line_number, line in enumerate(f):
            if line.startswith('>'):
                if line_number != 0:
                    log.warning('Second FASTA record detected - skipping.')
                    sys.exit(0)
            else:
                if line_number == 0:
                    log.error('FASTA does not begin with ">".')
                    sys.exit(1)
                else:
                    line = line.strip()
                    for base in line:
                        i += 1
                        if base != 'N':
                            bases.update(base)
                        if i == nbases:
                            get_bead(bases)
                            bases = Counter()
                            i = 0

        if i > 0:
            get_bead(bases)


def get_bead(bases):
    if sum(bases.values()) != 0:
        bead = bases.most_common(1)[0][0]
    else:
        bead = 'N'
    sys.stdout.write(f'{bead}\n')


if __name__ == '__main__':
    log = pct.create_logger()
    start = timer()
    RC = main()
    end = timer()
    log.info(f'Total time elapsed: {end - start} seconds.')
    sys.exit(RC)
