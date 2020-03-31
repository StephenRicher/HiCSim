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
from collections import defaultdict
from utilities import Atom, load_XYZ


def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__, infile=True,
        in_type='XYZ')
    parser.set_defaults(function=get_contact_frequency)

    parser.add_argument(
        '-x', '--xsize', required=True, type=float,
        help='Size of box in x dimension.')
    parser.add_argument(
        '-y', '--ysize', required=True, type=float,
        help='Size of box in y dimension.')
    parser.add_argument(
        '-z', '--zsize', required=True, type=float,
        help='Size of box in z dimension.')
    parser.add_argument(
        '--outdata', default='contacts.txt',
        help='Contact matrix output (default: %(default)s)')

    return (pct.execute(parser))


def get_contact_frequency(
        infile: str, outdata: str,
        xsize: float, ysize: float, zsize: float) -> None:

    log = pct.create_logger()
    contacts = defaultdict(lambda: defaultdict(int))
    xyz = load_XYZ(infile)
    for entry in xyz:
        log.info(f'Processing: {entry["comment"]}')
        contacts = update_contact_map(
            entry['atoms'], xsize, ysize, zsize, contacts)

    contacts = pd.concat(
        {k: pd.DataFrame.from_dict(v, 'index') for k, v in contacts.items()},
        axis=0, names=['atom1', 'atom2']).unstack(fill_value=0).to_numpy()

    np.savetxt(outdata, contacts)


def update_contact_map(
        atoms: List[Atom], xsize: float, ysize: float, zsize: float, contacts):

    for a, atom1 in enumerate(atoms):
        for b, atom2 in enumerate(atoms):
            if a == b:
                continue
            xdist = linear_distance(
                atom1.x, atom2.x, size=xsize, periodic=True)
            ydist = linear_distance(
                atom1.y, atom2.y, size=ysize, periodic=True)
            zdist = linear_distance(
                atom1.z, atom2.z, size=zsize, periodic=True)
            distance = distance_between(xdist, ydist, zdist)
            if distance < 3:
                contacts[a][b] += 1

    return contacts


def distance_between(xdist: float, ydist: float, zdist: float) -> float:
    return math.sqrt(xdist**2 + ydist**2 + zdist**2)


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
