#!/usr/bin/env python3

""" Mask a genome with relevant bases """

import os
import sys
import argparse
import pyCommonTools as pct
from utilities import commaPair
from contextlib import ExitStack
from subprocess import PIPE, Popen
from tempfile import NamedTemporaryFile

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__,)
    parser.set_defaults(function=mask)
    parser.add_argument(
        '--genome', metavar='FASTA', required=True,
        help='Input genome to perform masking.')
    parser.add_argument(
        '--bed', metavar='BED,CHAR', type=commaPair,
        action='append', required=True,
        help='BED file of regions to mask, paired with masking character.'
        'Call multiple times to add more files.')

    return (pct.execute(parser))


def mask(genome, bed):
    log = pct.create_logger()
    if duplicateMask(bed): sys.exit(1)

    allMasked = maskAllBases(genome)
    processes = []
    for i, (file, char) in enumerate(bed):
        input = allMasked if i == 0 else '/dev/stdin'
        cmd = ['bedtools', 'maskfasta', '-bed', file, '-mc', char,
               '-fi', input, '-fo', '/dev/stdout']
        stdin = None if i == 0 else processes[-1].stdout
        stdout = None if i == len(bed) - 1 else PIPE
        processes.append(Popen(cmd, stdin=stdin, stdout=stdout))

    for i, process in enumerate(processes):
        if i == len(processes) - 1:
            process.wait()
        else:
            process.stdout.close()
    os.remove(allMasked)


def duplicateMask(bed):
    """ Return True if 'bed' argument contains duplicates of masking character
        or the same file.
    """
    log = pct.create_logger()
    for index in [0, 1]:
        unique_entries = set([i[index] for i in bed])
        if len(unique_entries) != len(bed):
            if index == 0:
                log.error(f'{bed} contains duplicate filenames.')
            else:
                log.error(f'{bed} contains duplicate masking characters.')
            break
    else:
        return False
    return True


def maskAllBases(genome):
    """ Replace all bases in genome sequence with N. """

    with ExitStack() as stack:
        genome = stack.enter_context(open(genome))
        masked = stack.enter_context(
            NamedTemporaryFile(mode='wt', delete=False))
        for line in genome:
            if line.startswith('>'):
                masked.write(line)
            else:
                masked.write(f'{"N"*len(line.strip())}\n')
    return masked.name






if __name__ == '__main__':
    sys.exit(main())
