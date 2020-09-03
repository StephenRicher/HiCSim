#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency matrix """

import sys
import logging
import argparse
import numpy as np
from utilities import read_XYZ
from scipy.spatial.distance import cdist

__version__ = '1.0.0'

def main(dnaXYZ: str, monomerXYZ: str, out: str, distance: float, **kwargs) -> None:

    with open(dnaXYZ) as dna_fh, open(monomerXYZ) as monomer_fh:
        allDistances = []
        while True:
            try:
                dna = read_XYZ(dna_fh)
                monomer = read_XYZ(monomer_fh)
                # Compute distance of monomers to each DNA bead
                result = cdist(dna['atoms'], monomer['atoms'], 'euclidean')
                # Find closest monomer to each bead
                minDistance = np.amin(result, axis=1) < distance
                allDistances.append(minDistance)
            except EOFError:
                break
        # Convert list of distance arrays to 2D numpy matrix
        allDistances = np.vstack(tuple(allDistances))

    with open(out, 'wb') as f:
        np.save(f, allDistances)

def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'dnaXYZ', metavar='DNA',
        help='DNA coordinates in XYZ format')
    custom.add_argument(
        'monomerXYZ', metavar='MONOMER',
        help='Monomer coordinates in XYZ format')
    custom.add_argument(
        '--distance', default=1.3, type=float,
        help='Distance threshold for active transcription (default: %(default)s)')
    custom.add_argument(
        '--out', default='TFDistances.npy',
        help='Contact matrix output (default: %(default)s)')
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