#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency heatmap """

import sys
import pyCommonTools as pct
from typing import IO, List
from utilities import Atom, load_XYZ, print_XYZ, equal_n_atoms


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True,
        in_type='XYZ')
    parser.set_defaults(function=smooth_xyz)

    parser.add_argument(
        '-n', '--size', default=5, type=int,
        help='Number of timesteps to average.')
    parser.add_argument(
        '-o', '--overlap', default=0, type=int,
        help='Overlap between timestep averages.')

    return (pct.execute(parser))


def smooth_xyz(infile, size, overlap):

    log = pct.create_logger()
    if overlap >= size:
        log.error(
            f'Overlap ({overlap}) must be less than window size ({size}).')
        sys.exit(1)

    xyz = load_XYZ(infile)
    if not equal_n_atoms(xyz):
        log.error('Not all entries (timesteps) contain same number of atoms.')
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
        #print_mean_xyz(xyz, start, end)
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


if __name__ == '__main__':
    sys.exit(main())
