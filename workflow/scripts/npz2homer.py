#!/usr/bin/env python3

""" Write HiC matrix to Homer format """

import sys
import logging
import argparse
import pandas as pd
from utilities import npz
from scipy.sparse import load_npz


__version__ = '1.0.0'


def main(matrix, binsize, chromosome, start,  **kwargs) -> None:

    matrix = load_npz(matrix).toarray()
    names = [f'{chromosome}-{start+(i*binsize)}' for i in range(len(matrix))]
    matrix = pd.DataFrame(matrix, columns=names)
    matrix.insert(0, 'Regions', names)
    matrix.insert(0, 'HiCMatrix (directory=.)', names)
    matrix.to_csv(sys.stdout, sep='\t', index=False)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('matrix', type=npz,
        help='Contact matrix in .npz format.')
    requiredNamed = custom.add_argument_group(
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
    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')
    base.add_argument(
        '--verbose', action='store_const', const=logging.DEBUG,
        default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom])
    args = parser.parse_args()

    log_format='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    logging.basicConfig(level=args.verbose, format=log_format)

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
