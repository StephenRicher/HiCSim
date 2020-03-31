#!/usr/bin/env python3

import pyCommonTools as pct
from timeit import default_timer as timer


class Atom:

    def __init__(self, record=None):
        if record is None:
            self._record = ['.', 0, 0, 0]
        else:
            self._record = record.strip().split()

    def __repr__(self):
        return ' '.join([str(i) for i in self._record])

    @property
    def type(self):
        return self._record[0]

    @type.setter
    def type(self, var):
        self._record[0] = var

    @property
    def x(self):
        return float(self._record[1])

    @x.setter
    def x(self, var):
        self._record[1] = float(var)

    @property
    def y(self):
        return float(self._record[2])

    @y.setter
    def y(self, var):
        self._record[2] = float(var)

    @property
    def z(self):
        return float(self._record[3])

    @z.setter
    def z(self, var):
        self._record[3] = float(var)


def load_XYZ(XYZ_path: str):
    """ Read XYZ file and parse into list of dictionaries. """

    log = pct.create_logger()
    log.debug(f'Reading {XYZ_path} into memory.')
    start = timer()
    xyz = []
    with open(XYZ_path) as f:
        for id, line in enumerate(f):
            n_atoms = int(line.strip())
            comment = next(f).strip()
            xyz.append(
                {'n_atoms' : n_atoms,
                 'comment' : comment,
                 'atoms'   : []})
            for n in range(n_atoms):
                xyz[id]['atoms'].append(Atom(next(f)))
    end = timer()
    log.debug(f'Parsed {XYZ_path} in {end - start} seconds.')

    return xyz


def print_XYZ(xyz) -> None:
    """ Write XYZ object to stdout. """

    for entry in xyz:
        print(entry['n_atoms'])
        print(entry['comment'])
        for atom in entry['atoms']:
            print(atom)


def equal_n_atoms(xyz):
    """ Return TRUE if same number of atoms across all entries (timesteps) """

    return all([entry['n_atoms'] for entry in xyz])
