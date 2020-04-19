#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency matrix """

import sys
import numpy as np
import pyCommonTools as pct
from utilities import read_XYZ
from timeit import default_timer as timer
from scipy.spatial.distance import pdist, squareform

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
        '--outdata', default='contacts.txt',
        help='Contact matrix output (default: %(default)s)')

    return (pct.execute(parser))


def get_contact_frequency(infile: str, outdata: str, distance: float) -> None:

    log = pct.create_logger()
    start = timer()

    sqdistance = distance**2
    contacts = 0
    with open(infile) as f:
        while True:
            try:
                xyz = read_XYZ(f)
                contacts += pdist(xyz['atoms'], 'sqeuclidean') < sqdistance
                end = timer()
            except EOFError:
                break

    np.savetxt(outdata, squareform(contacts))

    end = timer()
    log.info(f'Contact matrix created in {end - start} seconds.')


if __name__ == '__main__':
    sys.exit(main())
