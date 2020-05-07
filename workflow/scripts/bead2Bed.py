#!/usr/bin/env python3

""" Read bead file and output to BED9 with RGB colours for beads """

import sys
import pyCommonTools as pct


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__,
        infile=True, in_type='BEADS')
    parser.set_defaults(function=bead2Bed)

    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--binsize', type=int, required=True,
        help='Number of bases used to represent each bin.')
    requiredNamed.add_argument(
        '--chromosome', required=True,
        help='Chromosome of matrix.')
    requiredNamed.add_argument(
        '--start', type=int, required=True,
        help='Start coordinate of matrix.')

    return (pct.execute(parser))


def bead2Bed(infile, binsize, chromosome, start) -> None:

    rgb = {'F' : '255,0,0',
           'R' : '0,0,255'}

    position = start - 1
    with open(infile) as fh:
        for line in fh:
            bead = line.strip()
            if bead != 'N':
                start = position
                end = position + binsize
                name = '.'
                score = 1000
                strand = '.'
                thickStart = start
                thickEnd = start
                try:
                    itemRgb = rgb[bead]
                except KeyError:
                    itemRgb = '0,0,0'
                print(chromosome, start, end, name, score,
                    strand, thickStart, thickEnd, itemRgb, sep='\t')
            position = position + binsize


if __name__ == '__main__':
    sys.exit(main())
