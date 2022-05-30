#!/usr/bin/env python3

""" Write HiC matrix to Homer format """


import sys
import argparse
import pandas as pd
from scipy import ndimage
from scipy.sparse import load_npz
from hicmatrix import HiCMatrix as hm
from utilities import setDefaults, createMainParent


__version__ = '1.0.0'


def npz2homer(npz: str, expHiC: str, binSize: int, chromosome: str, start: int) -> None:

    matrix = load_npz(npz).toarray()
    if expHiC is not None:
        hic = hm.hiCMatrix(expHiC)
        zoom = binSize / hic.getBinSize()
        matrix = ndimage.zoom(matrix, zoom).clip(min=0)
        depthFactor = hic.matrix.sum() / matrix.sum()
        matrix = (matrix * depthFactor).astype(int)
        binSize = hic.getBinSize()

    names = [f'{chromosome}-{start+(i*binSize)}' for i in range(len(matrix))]
    matrix = pd.DataFrame(matrix, columns=names)
    matrix.insert(0, 'Regions', names)
    matrix.insert(0, 'HiCMatrix (directory=.)', names)
    # Set NaN (diagonal) to zero and round scientific notation numbers
    matrix = matrix.fillna(0).round(10)
    matrix.to_csv(sys.stdout, sep='\t', index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=npz2homer)
    parser.add_argument('npz', metavar='NPZ', help='Input numpy matrix.')
    parser.add_argument(
        '--expHiC',
        help='Experimental HiC matrix in .h5/.cool '
             'format to resize simulation matrix against.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--binSize', type=int, required=True,
        help='Number of bases to represent each bin.')
    requiredNamed.add_argument(
        '--chromosome', required=True, help='Chromosome of matrix.')
    requiredNamed.add_argument(
        '--start', type=int, required=True, help='Start coordinate of matrix.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
