#!/usr/bin/env python3

import sys
import pyCommonTools as pct
from contextlib import ExitStack

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True, in_type='BED')
    parser.add_argument(
        '-r', '--min_rep', type=int, default=1,
        help='Minimum number of replicates required per record to write.')
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--reverse', required=True,
        help='Output file for reverse orientation intervals.')
    requiredNamed.add_argument(
        '--forward', required=True,
        help='Output file for forward orientation intervals.')

    parser.set_defaults(function=splitbed)

    return (pct.execute(parser))


def splitbed(infile, forward, reverse, min_rep):

    with ExitStack() as stack:

        f = stack.enter_context(pct.open(infile))
        out_forward = stack.enter_context(pct.open(forward, 'w'))
        out_reverse = stack.enter_context(pct.open(reverse, 'w'))

        for line in f:
            columns = line.strip().split()
            rep = int(columns[3])
            if rep < min_rep:
                continue
            orientation = columns[5]
            if orientation == '+':
                out = out_forward
            else:
                out = out_reverse

            print(line, end='', file = out)


if __name__ == "__main__":
    main()
