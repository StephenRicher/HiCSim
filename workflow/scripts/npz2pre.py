#!/usr/bin/env python3

""" Write HiC matrix to short (score) format for Juicer Pre """

import sys
import pandas as pd
from utilities import npz
import pyCommonTools as pct
from scipy.sparse import load_npz

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(verbose=True, version=__version__)
    parser.set_defaults(function=npz2pre)

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


def npz2pre(matrix, binsize, chromosome, start) -> None:

    matrix = load_npz(matrix).toarray()
    names = [f'{start+(i*binsize)}' for i in range(len(matrix))]
    matrix = pd.DataFrame(matrix, index=names, columns=names)
    matrix = matrix.stack().reset_index()
    matrix.insert(0, 'str1', 0)
    matrix.insert(1, 'chr1', chromosome)
    matrix.insert(3, 'frag1', 0)
    matrix.insert(4, 'str2', 1)
    matrix.insert(5, 'chr2', chromosome)
    matrix.insert(7, 'frag2', 1)
    matrix.to_csv(sys.stdout, sep='\t', index=False, header=False)


if __name__ == '__main__':
    sys.exit(main())
