#!/usr/bin/env python3

""" Run CTCF motif prediction API """

import os
import sys
import logging
import requests
import argparse
import fileinput
import pandas as pd

__version__ = '1.0.0'


def main(file, **kwargs):
    """ CTCF prediction """

    sequence = ''
    with fileinput.input(file) as fh:
        for line in fh:
            if sequence and line.startswith('>'):
                predictCTCF(sequence)
                sequence = line
            else:
                sequence += line
        predictCTCF(sequence)


def predictCTCF(sequence):
    data= {'my_seqData': sequence, 'my_run' : 1}
    r = requests.post('http://insulatordb.uthsc.edu/storm_new.php', data=data)
    results_table = pd.read_html(r.text)[2]
    results_tsv = results_table.to_csv(index=False, header=False, sep='\t')
    if len(results_tsv) == '':
        sys.stderr.write(sequence)
    else:
        sys.stdout.write(results_tsv)


def parse_arguments(
        version=None, verbose=True,
        epilog='Stephen Richer, University of Bath, Bath, UK '
               '(sr467@bath.ac.uk)',
        log_format='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'):

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='FASTA', nargs='?', default=[],
        help='Fasta file of putative CTCF sites (default: stdin)')

    base = argparse.ArgumentParser(add_help=False)
    if version:
        base.add_argument(
            '--version', action='version', version=f'%(prog)s {version}')
    if verbose:
        base.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom])
    args = parser.parse_args()

    verbosity = args.verbose if verbose else None
    logging.basicConfig(level=verbosity, format=log_format)

    return args


if __name__ == '__main__':
    args = parse_arguments(version=__version__)
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
