#!/usr/bin/env python3

""" Generate synthetic polymer """

import sys
import math
import random
import logging
import argparse

__version__ = '1.0.0'


def main(subclusterSequence, nClusters, bead, buffer, intraClusterBeads,
        interClusterDistance, intraClusterDistance, **kwargs):

    writeBeads(buffer)
    for cluster in range(1, nClusters+1):
        writeSubCluster(subclusterSequence)
        writeIntraClusterBeads(intraClusterDistance, intraClusterBeads, bead)
        writeSubCluster(subclusterSequence)
        # If last cluster don't write interClusterDistance
        if cluster < nClusters:
            writeBeads(interClusterDistance, bead)
    writeBeads(buffer)


def writeIntraClusterBeads(intraClusterDistance, intraClusterBeads, bead='N'):
    """ Distribute TU beads with even gap distance """

    # Sequence begins and ends with gap so number of gaps is 1 plus TU beads
    nGaps = intraClusterBeads + 1
    gapLength = math.floor((intraClusterDistance - intraClusterBeads) / nGaps)
    # Extra gap length remaining
    extraLength = (intraClusterDistance - intraClusterBeads) % nGaps
    # Assign extra gap length to central gap
    # If gaps is even - randomly assign central gap
    if (nGaps % 2) == 0:
        center = nGaps / 2
        middleGap = random.choice([center - 1, center + 1])
    else:
        middleGap = math.ceil(nGaps / 2)
    for gap in range(1, nGaps + 1):
        if gap == middleGap:
            nGapBeads = gapLength + extraLength
        else:
            nGapBeads = gapLength
        writeBeads(nGapBeads)
        if gap < nGaps:
            print('3')


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
        'subclusterSequence', nargs='?', default='3N3N3N3',
        help='Sequence of sub-cluster')
    custom.add_argument(
        '--nClusters',type=int, default=3,
        help='Number of clusters in polymer')
    custom.add_argument(
        '--interClusterDistance', type=int, default=349,
        help='Distance between clusters')
    custom.add_argument(
        '--intraClusterDistance', type=int, default=70,
        help='Distance between subclusters')
    custom.add_argument(
        '--intraClusterBeads', type=int, default=4,
        help='Number of beads between the sub-clusters.')
    custom.add_argument(
        '--buffer', type=int, default=25,
        help='Number of edge beads to polymer.')

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
