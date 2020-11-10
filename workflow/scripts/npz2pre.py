#!/usr/bin/env python3

""" Write HiC matrix to short (score) format for Juicer Pre """


import sys
import pandas as pd
from utilities import setDefaults
from scipy.sparse import load_npz


__version__ = '1.0.0'


def npz2pre(npz: str, binsize: int, chromosome: str, start: int) -> None:

    matrix = load_npz(npz).toarray()
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


def parseArgs():

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('npz', metavar='NPZ', help='Input numpy matrix.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binsize', type=int, required=True,
        help='Number of bases used to represent each bin.')
    requiredNamed.add_argument(
        '--chromosome', required=True,
        help='Chromosome of matrix.')
    requiredNamed.add_argument(
        '--start', type=int, required=True,
        help='Start coordinate of matrix.')

    return setDefaults(parser, verbose=False, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(npz2pre(**vars(args)))
