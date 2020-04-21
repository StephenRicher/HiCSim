#!/usr/bin/env python3

import re
import argparse
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


def read_XYZ(XYZ_fobj):
    """ Read a timepoint of XYZ coordinates each time called. """

    line = XYZ_fobj.readline()
    if not line:
        raise EOFError
    n_atoms = int(line.strip())
    comment = next(XYZ_fobj).strip()
    xyz = {'n_atoms' : n_atoms,
           'comment' : comment,
           'atoms'   : []}
    for n in range(n_atoms):
        atom = Atom(next(XYZ_fobj))
        xyz['atoms'].append([atom.x, atom.y, atom.z])

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


def npz(value):
    ''' Ensure output file ends with '.npz'. '''

    if not value.endswith('.npz'):
        raise argparse.ArgumentTypeError(
            f'Out file {value} must end with: \'.npz\'.')
    else:
        return value


def region(value):
    ''' Validate input for region argument of create_contact_matrix '''

    pattern = re.compile('^[0-9]+-[0-9]+$')
    if not pattern.match(value):
        raise argparse.ArgumentTypeError(
            f'Expected format is START-END e.g 1-100000000')
    regions = {}
    regions['start'] = int(value.split('-')[0])
    regions['end'] = int(value.split('-')[1])

    if not regions['start'] < regions['end']:
        raise argparse.ArgumentTypeError(
            f'Start coordinate {regions["start"]} not less '
            f'than end coordinate {regions["end"]}.')
    else:
        return regions
