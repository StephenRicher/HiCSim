#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency matrix """

import sys
import math
import numpy as np
import pandas as pd
import seaborn as sns
import pyCommonTools as pct
from typing import IO, List
import matplotlib.pyplot as plt
from itertools import combinations, repeat
from utilities import Atom, load_XYZ
from timeit import default_timer as timer
import ctypes
import multiprocessing


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True,
        in_type='XYZ')
    parser.set_defaults(function=get_contact_frequency)

    parser.add_argument(
        '-d', '--distance', default=3, type=float,
        help='Max contact distance between particles (default: %(default)s)')
    parser.add_argument(
        '-x', '--xsize', required=True, type=float,
        help='Size of box in x dimension')
    parser.add_argument(
        '-y', '--ysize', required=True, type=float,
        help='Size of box in y dimension')
    parser.add_argument(
        '-z', '--zsize', required=True, type=float,
        help='Size of box in z dimension')
    parser.add_argument(
        '-@', '--threads', default=1, type=int,
        help='Number of threads for multiprocessing (default: %(default)s)')
    parser.add_argument(
        '--outdata', default='contacts.txt',
        help='Contact matrix output (default: %(default)s)')

    return (pct.execute(parser))


def get_contact_frequency(
        infile: str, outdata: str, distance: float, threads: int,
        xsize: float, ysize: float, zsize: float) -> None:

    log = pct.create_logger()
    start = timer()

    xyz = load_XYZ(infile)
    n_atoms = max(i['n_atoms'] for i in xyz)

    global contacts
    shared_array_base = multiprocessing.Array(ctypes.c_double, n_atoms*n_atoms)
    contacts = np.ctypeslib.as_array(shared_array_base.get_obj())
    contacts = contacts.reshape(n_atoms, n_atoms)

    max_square_distance = distance**2

    # Generate tuple of args for submission to update_contact_map()
    args = zip(xyz,repeat(xsize), repeat(xsize),
        repeat(xsize), repeat(max_square_distance))
    pool = multiprocessing.Pool(processes=threads)
    pool.starmap(update_contact_map, args)

    # Copy upper triangle to lower triangle
    contacts = contacts + contacts.T - np.diag(np.diag(contacts))
    np.savetxt(outdata, contacts)

    end = timer()
    log.info(f'Contact matrix created in {end - start} seconds.')


# Also run in parallel?
def update_contact_map(
        atoms: List[Atom], xsize: float, ysize: float, zsize: float,
        max_square_distance: float):

    # Loop through unique atom pairs - don't run same pairs twice!
    for A1, A2 in combinations(enumerate(atoms['atoms']), 2):

        atom1 = A1[1]
        atom2 = A2[1]
        xdist = linear_distance(
            atom1.x, atom2.x, size=xsize, periodic=True)
        ydist = linear_distance(
            atom1.y, atom2.y, size=ysize, periodic=True)
        zdist = linear_distance(
            atom1.z, atom2.z, size=zsize, periodic=True)
        square_distance = xdist**2 + ydist**2 + zdist**2
        if square_distance < max_square_distance:
            contacts[A1[0]][A2[0]] += 1

    return contacts


def linear_distance(
        p1: float, p2: float, size: float, periodic: bool = False) -> float:
    """ Calculate distance between two atoms in 1 dimension. """

    distance = abs(p1 - p2)
    if periodic:
        return min(size - distance, distance)
    else:
        return distance


if __name__ == '__main__':
    sys.exit(main())
