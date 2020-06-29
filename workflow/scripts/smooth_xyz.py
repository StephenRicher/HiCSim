#!/usr/bin/env python3

""" Smooth XYZ atom positions across timepoints. """

import sys
import logging
import argparse
from typing import IO, List
from utilities import Atom, load_XYZ, print_XYZ, equal_n_atoms

__version__ = '1.0.0'

def main(file :str , size : int , overlap : int, **kwargs):

    if overlap >= size:
        log.error(
            f'Overlap ({overlap}) must be less than window size ({size}).')
        sys.exit(1)

    xyz = load_XYZ(file)
    if not equal_n_atoms(xyz):
        logging.error(
            'Not all entries (timesteps) contain same number of atoms.')
        sys.exit(1)

    n_timesteps = len(xyz)
    increment = size - overlap
    start_indexes = range(0, n_timesteps + increment, increment)
    end_indexes = range(size - 1, n_timesteps + increment, increment)
    final = False  # Ensure loop exits once end_index reaches max index

    mean_xyz = []
    for start, end in zip(start_indexes, end_indexes):
        if end >= n_timesteps - 1:
            final = True
            end = n_timesteps - 1
        mean_xyz.append(mean_xyz_entry(xyz, start, end))
        if final:
            break
    print_XYZ(mean_xyz)


def mean_xyz_entry(xyz, start, end):
    n_atoms = xyz[0]['n_atoms']
    comment = ''
    atoms = []
    for atom_id in range(n_atoms):
        mean_atom = Atom()
        size = end - start + 1
        for step in range(start, end + 1):
            atom = xyz[step]['atoms'][atom_id]
            if step == start:
                mean_atom.type = atom.type
            mean_atom.x += atom.x
            mean_atom.y += atom.y
            mean_atom.z += atom.z
        mean_atom.x *= (1/size)
        mean_atom.y *= (1/size)
        mean_atom.z *= (1/size)
        atoms.append(mean_atom)
    entry = (
        {'n_atoms' : n_atoms,
        'comment' : comment,
        'atoms'   : atoms})
    return entry


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'file', metavar='XYZ', nargs='?', default=[],
        help='Input XYZ file (default: stdin)')
    custom.add_argument(
        '--size', default=5, type=int,
        help='Number of timesteps to average (default: %(default)s)')
    custom.add_argument(
        '--overlap', default=0, type=int,
        help='Overlap between timestep averages (default: %(default)s)')
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
