#!/usr/bin/env python3

""" Generate compute N-masked FASTA sequence """

import os
import sys
import pyCommonTools as pct

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(version=__version__,)
    parser.set_defaults(function=maskAllBases)
    parser.add_argument(
        'fasta', metavar='FASTA', help='Input FASTA file to mask.')

    return (pct.execute(parser))


def maskAllBases(fasta):
    """ Replace all bases in genome sequence with N. """

    with pct.open(fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                sys.stdout.write(line)
            else:
                sys.stdout.write(f'{"N"*len(line.strip())}\n')


if __name__ == '__main__':
    sys.exit(main())
