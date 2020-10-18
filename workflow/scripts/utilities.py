#!/usr/bin/env python3

import sys
import re
import argparse
import pandas as pd
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

    return xyz


def read_XYZ(XYZ_fobj, includeID=[], excludeTypes=[]):
    """ Read a timepoint of XYZ coordinates each time called. """

    line = XYZ_fobj.readline()
    if not line:
        raise EOFError
    n_atoms = int(line.strip())
    comment = next(XYZ_fobj).strip()
    xyz = {'n_atoms' : n_atoms,
           'comment' : comment,
           'atoms'   : []}
    excluded = 0
    for n in range(1, n_atoms + 1):
        type, x, y , z = next(XYZ_fobj).strip().split()
        if (n not in includeID) or (type in excludeTypes):
            excluded += 1
        else:
            xyz['atoms'].append([x, y, z])
    xyz['n_atoms'] -= excluded

    return xyz


def readCustom(fobj, includeIDs=[], excludeTypes=[]):
    """ Read a timepoint of XYZ coordinates each time called. """

    # Convert user input to str to ensure match between file
    includeIDs = [str(id_) for id_ in includeIDs]
    excludeTypes = [str(type_) for type_ in excludeTypes]

    try:
        header = [next(fobj).strip() for line in range(9)]
    except StopIteration:
        raise EOFError
    xyz = {'timestep': int(header[1]),
           'nAtoms'  : int(header[3]),
           'atoms'   : []}
    excluded = 0
    for n in range(xyz['nAtoms']):
        id_, type_, x, y, z, ix, iy, iz = next(fobj).strip().split()
        if (includeIDs and id_ not in includeIDs) or (type_ in excludeTypes):
            excluded += 1
        else:
            xyz['atoms'].append([float(x), float(y), float(z)])
    xyz['nAtoms'] -= excluded

    return xyz


def readCustom2(fobj, includeIDs=[], excludeTypes=[]):
    """ Read a timepoint of XYZ coordinates into pandas. """

    try:
        header = [next(fobj).strip() for line in range(9)]
    except StopIteration:
        raise EOFError

    nAtoms = int(header[3])
    time = int(header[1])
    atomId, atomType, xPos, yPos, zPos = [], [], [], [], []
    for n in range(nAtoms):
        id_, type_, x, y, z, ix, iy, iz = next(fobj).strip().split()
        if (includeIDs and id_ not in includeIDs) or (type_ in excludeTypes):
            continue
        atomId.append(int(id_))
        atomType.append(type_)
        xPos.append(float(x))
        yPos.append(float(y))
        zPos.append(float(z))

    data = pd.DataFrame(
        {'id'   : atomId  ,
         'type' : atomType,
         'time' : time    ,
         'xPos' : xPos    ,
         'yPos' : yPos    ,
         'zPos' : zPos    })
    return data



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


def scale(value):
    ''' Ensure input is numeric or 'None'. '''

    if value == 'None':
        return None
    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f'must be float value or None')


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

def coordinates(value):
    ''' Validate input for genomic coordinates  '''

    pattern = re.compile('^[^:-]+:[0-9]+-[0-9]+$')
    if not pattern.match(value):
        raise argparse.ArgumentTypeError(
            'Expected format is CHR:START-END e.g chr1:1-1000. '
            'Chromosome name cannot contain ": -" .')
    coords = {}
    coords['chr'], coords['start'], coords['end'] = re.sub('[:-]', ' ', value).split()
    coords['start'] = int(coords['start'])
    coords['end'] = int(coords['end'])
    if not coords['start'] < coords['end']:
        raise argparse.ArgumentTypeError(
            f'Start coordinate {coords["start"]} not less '
            f'than end coordinate {coords["end"]}.')
    else:
        return coords


def commaPair(value):
    ''' Split command seperated pair and return as dictionary '''

    value = value.split(',')
    v1 = value[0].strip()
    v2 = value[1].strip()
    if len(v2) != 1:
        raise argparse.ArgumentTypeError(
            f'Masking character {v2} shoud be single character.')
    else:
        return (v1, v2)


def coeff(value):
    ivalue = float(value)
    if not -1 <= ivalue <= 1:
        raise argparse.ArgumentTypeError(f'{value} must be between -1 and 1.')
    return ivalue
