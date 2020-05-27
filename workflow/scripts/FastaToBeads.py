#!/usr/bin/env python3

""" Script to read a single record FASTA file and compress N bases to 1 bead """

import sys
import random
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(version=__version__)
    parser.add_argument('fasta', metavar='FASTA', help='Input FASTA file')
    parser.add_argument(
        '--nbases', default=1000, type=int,
        help='Number of base to represent 1 bead.')
    parser.set_defaults(function=compress)

    return (pct.execute(parser))


def compress(fasta, nbases):

    log = pct.create_logger()
    with pct.open(fasta) as fh:
        bases = []
        for i, line in enumerate(fh):
            if line.startswith('>'):
                if i != 0:
                    log.warning('Second FASTA record detected - skipping.')
                    break
            else:
                if i == 0:
                    log.error('FASTA does not begin with ">".')
                    sys.exit(1)
                else:
                    for base in line.strip():
                        bases.append(base)
                        if len(bases) == nbases:
                            sys.stdout.write(f'{get_bead(bases)}\n')
                            bases = []
        if len(bases) > 0:
            sys.stdout.write(f'{get_bead(bases)}\n')


def get_bead(bases):
    # Remove any 'N' bases from list
    bases = [base for base in bases if base != 'N']
    if len(bases) == 0:
        bead = 'N'
    else:
        # If both CTCF orientations exist create a bead 'B' to represent both
        if 'F' in bases and 'R' in bases:
            bases = ['B' if base in ['F', 'R'] else base for base in bases]
        bead = random.choice(bases)
    return bead


if __name__ == '__main__':
    sys.exit(main())
