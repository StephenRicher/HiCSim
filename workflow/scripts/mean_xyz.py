#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency heatmap """

import sys
import pyCommonTools as pct
from typing import IO, List
from timeit import default_timer as timer
from utilities import Atom, load_XYZ, print_XYZ


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__)
    parser.set_defaults(function=compute_mean_xyz)

    parser.add_argument(
        'XYZs', nargs='+',
        help='Input XYZ files to merge')

    return (pct.execute(parser))


def compute_mean_xyz(XYZs):

    log = pct.create_logger()
    log.info(f'Begining averaging for XYZ files.')
    start = timer()
    n_files = len(XYZs)
    for i, xyz in enumerate(XYZs):
        if i == 0:
            summed_xyz = load_XYZ(xyz)
        else:
            summed_xyz = sum_XYZ(summed_xyz, load_XYZ(xyz))
    mean_xyz = scale_XYZ(summed_xyz, 1/n_files)
    print_XYZ(mean_xyz)
    end = timer()
    log.info(f'Averaging completed in {end - start} seconds.')

def sum_XYZ(xyz1, xyz2):
    """ Sum the X, Y and Z coordinate of atoms in a pair of XYZ objects. """


    log = pct.create_logger()
    log.debug(f'Summing XYZ files.')
    start = timer()
    if len(xyz1) != len(xyz2):
        log.error('Number of entries between XYZ files does not match.')
        sys.exit(1)

    summed_xyz = []
    for entry1, entry2 in zip(xyz1, xyz2):

        if entry1['comment'] != entry2['comment']:
            log.warning(
                f'XYZ Comments {entry1["comment"]} and {entry2["comment"]} '
                'do not match for same entry. Using first comment.')
        if entry1['n_atoms'] != entry2['n_atoms']:
            log.error('Number of atoms between XYZ files does not match.')
            sys.exit(1)

        atoms = []
        for atom1, atom2 in zip(entry1['atoms'], entry2['atoms']):
            if atom1.type != atom2.type:
                log.error('Atom type between XYZ files does not match.')
                sys.exit(1)
            sum_atom = Atom()
            sum_atom.type = atom1.type
            sum_atom.x = atom1.x + atom2.x
            sum_atom.y = atom1.y + atom2.y
            sum_atom.z = atom1.z + atom2.z
            atoms.append(sum_atom)

        entry = (
            {'n_atoms' : entry1['n_atoms'],
             'comment' : entry1['comment'],
             'atoms'   : atoms})
        summed_xyz.append(entry)
    end = timer()
    log.debug(f'Summed XYZ files in {end - start} seconds.')

    return summed_xyz


def scale_XYZ(xyz, scale_factor: int = 1):
    """ Scale XYZ coordinates """


    log = pct.create_logger()
    log.info(f'Rescaling XYZ by a scale factor of {scale_factor}.')
    start = timer()
    for entry in xyz:
        for atom in entry['atoms']:
            atom.x *= scale_factor
            atom.y *= scale_factor
            atom.z *= scale_factor
    end = timer()
    log.info(f'Rescaled in {end - start} seconds.')

    return xyz



if __name__ == '__main__':
    sys.exit(main())
