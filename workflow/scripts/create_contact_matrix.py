#!/usr/bin/env python3

""" Script to read LAMMPS output and generate contact frequency matrix """

import sys
import numpy as np
import pyCommonTools as pct
from utilities import read_XYZ, npz
from timeit import default_timer as timer
from scipy.sparse import save_npz, csc_matrix
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
        '--outdata', default='contacts.npz', type=npz,
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
                start_step = timer()
                xyz = read_XYZ(f)
                contacts += pdist(xyz['atoms'], 'sqeuclidean') < sqdistance
                end_step = timer()
                log.info(f'Block processed in {end_step - start_step} seconds.')
            except EOFError:
                break

    save_npz(outdata, csc_matrix(squareform(contacts)))

    end = timer()
    log.info(f'Contact matrix created in {end - start} seconds.')


if __name__ == '__main__':
    sys.exit(main())
