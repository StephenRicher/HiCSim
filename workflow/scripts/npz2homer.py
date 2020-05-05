#!/usr/bin/env python3

""" Write HiC matrix to Homer format """

import sys
import pandas as pd
from utilities import npz
import pyCommonTools as pct
from scipy.sparse import load_npz

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)
    parser.set_defaults(function=npz2homer)

    parser.add_argument('matrix', type=npz)
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--binsize', type=int, required=True,
        help='Number of bases used to represent each bin.')
    requiredNamed.add_argument(
        '--chromosome', required=True,
        help='Chromosome of matrix.')
    requiredNamed.add_argument(
        '--start', type=int, required=True,
        help='Start coordinate of matrix.')

    return (pct.execute(parser))


def npz2homer(matrix, binsize, chromosome, start) -> None:

    matrix = load_npz(matrix).toarray()
    names = [f'{chromosome}-{start+(i*binsize)}' for i in range(len(matrix))]
    matrix = pd.DataFrame(matrix, columns=names)
    matrix.insert(0, 'Regions', names)
    matrix.insert(0, 'HiCMatrix (directory=.)', names)
    matrix.to_csv(sys.stdout, sep='\t', index=False)


if __name__ == '__main__':
    sys.exit(main())
