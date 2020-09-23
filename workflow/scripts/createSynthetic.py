#!/usr/bin/env python3

""" Generate synthetic polymer """

import sys
import logging
import argparse

__version__ = '1.0.0'


def main(subclusterSequence, nClusters, bead,
        interClusterDistance, subclusterDistance, **kwargs):

    for cluster in range(nClusters):
        writeSubCluster(subclusterSequence)
        writeBeads(subclusterDistance, bead)
        writeSubCluster(subclusterSequence)
        if cluster != nClusters - 1:
            writeBeads(interClusterDistance, bead)


def writeSubCluster(sequence):
    for bead in sequence:
        print(bead)


def writeBeads(nBeads, bead='N'):
    for i in range(nBeads):
        print(bead)

def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'subclusterSequence', nargs='?', default='3NNNN3NNNN3NNN3NNNN3NNNN3',
        help='Sequence of sub-cluster')
    custom.add_argument(
        '--nClusters',type=int, default=3,
        help='Number of clusters in polymer')
    custom.add_argument(
        '--interClusterDistance', type=int, default=350,
        help='Distance between clusters')
    custom.add_argument(
        '--subclusterDistance', type=int, default=50,
        help='Distance between subclusters')
    custom.add_argument(
        '--bead', default='N',
        help='Default bead for interval sequences.')
    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')
    base.add_argument(
        '--verbose', action='store_const', const=logging.DEBUG,
        default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom])
    args = parser.parse_args()

    log_format='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    logging.basicConfig(level=args.verbose, format=log_format)

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
