#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency heatmap """

import sys
import math
import pyCommonTools as pct
from typing import IO, List
from collections import defaultdict
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

    return (pct.execute(parser))


class Atom:

    def __init__(self, record):
        self.record = record.strip().split()

    def __repr__(self):
        return ' '.join(self.record)

    @property
    def type(self):
        return self.record[0]

    @property
    def x(self):
        return float(self.record[1])

    @property
    def y(self):
        return float(self.record[2])

    @property
    def z(self):
        return float(self.record[3])


def get_contact_frequency(
        infile: str, xsize: float, ysize: float, zsize: float):

    contacts = defaultdict(lambda: defaultdict(int))
    with pct.open(infile) as f:
        n_atoms = int(next(f))
        for line in f:
            if 'Timestep' in line:
                sys.stderr.write(f'Processing: {line}')
                atoms = load_timestep(f, n_atoms)
                contacts = update_contact_map(atoms, xsize, ysize, zsize, contacts)

    contacts = pd.concat({k: pd.DataFrame.from_dict(v, 'index') for k, v in contacts.items()},
        axis=0, names=['atom1', 'atom2']).unstack(fill_value=0).to_numpy()
    np.savetxt('/home/stephen/Desktop/contacts.txt', contacts)
    sns.heatmap(np.log10(contacts + 1))
    #plt.show()
    plt.savefig('/home/stephen/Desktop/contacts.png', dpi=300)



def load_timestep(f: IO, n_atoms: int) -> List[Atom]:
    """ Return a list of Atoms for a timestep. """

    atoms = []
    for atom in range(1, n_atoms + 1):
        try:
            atoms.append(Atom(next(f)))
        except StopIteration:
            # Add warning here
            break
    return atoms


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
