#!/usr/bin/env python3

""" Randomly insert beads of specific type into bead sequence  """


import re
import sys
import argparse
import pandas as pd
from argUtils import setDefaults, createMainParent, positiveInt


__version__ = '1.0.0'


def addTUbeads(nBeads: int, sequence: str, bead: str, seed: float):
    beadSequence = pd.read_csv(sequence, header=None)[0]
    # Get index positions of all N beads
    Nbases = pd.Series(beadSequence[beadSequence == 'N'].index)
    # Randomly select beads to convert to TU beads
    try:
        TUbeads = Nbases.sample(nBeads, random_state=seed)
    except ValueError:
        logging.error(
            f'Number of TUs to add ({nBeads}) is higher than '
            f'number of N beads ({lenNbases}) in sequence.')
    beadSequence[TUbeads] = bead
    # Write to stdout
    beadSequence.to_csv(sys.stdout, header=False, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=addTUbeads)
    parser.add_argument(
        'nBeads', type=positiveInt,
        help='Number of TU beads to insert in sequence.')
    parser.add_argument(
        'sequence', nargs='?', default=sys.stdin,
        help='Bead sequence to add (default: stdin)')
    parser.add_argument(
        '--bead', default='3', type=bead,
        help='Single character to use as bead (default: %(default)s)')
    parser.add_argument(
        '--seed', default=None, type=float,
        help='Seed for determing TU positions (default: %(default)s)')

    return setDefaults(parser)


def bead(value):
    invalid = 'NRFB'
    if not re.match(rf'^[^{invalid}]$', value):
        raise argparse.ArgumentTypeError(
            f'{value} must be single character and cannot be any of {invalid}.')
    return value


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
