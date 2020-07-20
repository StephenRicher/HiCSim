#!/usr/bin/env python3

""" Plot Radius of Gyration across replicates from LAMMPS output """

import sys
import logging
import argparse
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

__version__ = '1.0.0'

def main(files, dpi, confidence, format, output, **kwargs):

    data = np.dstack([np.loadtxt(rep) for rep in files])
    mean = np.mean(data, axis=2)
    time = mean[:,0]
    mean = mean[:,1]
    n = len(files)
    std_err = stats.sem(data, axis=2)[:,1]
    ci = std_err * stats.t.ppf((1 + confidence) / 2, n - 1)
    dydx = np.diff(mean)/ np.diff(time)
    midpoints = (time[1:] + time[:-1]) / 2

    fig, (ax1, ax2) = plt.subplots(2,1)

    ax1.plot(time, mean)
    ax1.fill_between(time, (mean - ci), (mean + ci), color='b', alpha=.1)
    ax1.axvline(x=1, color='Red', alpha=0.3)
    ax1.tick_params(labelsize=16)
    ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel('RG', fontsize=20)
    ax1.set_title(
        f'Radius of Gyration over time with {confidence:.0%} CI (n = {n})',
        loc='left', fontsize=20)

    ax2.plot(midpoints, dydx)
    ax2.axvline(x=1, color='Red', alpha=0.3)
    ax2.axhline(y=0, color='Red', alpha=0.3)
    ax2.tick_params(labelsize=16)
    ax2.set_xlabel('time', fontsize=20)
    ax2.set_ylabel('dRG/dt', fontsize=20)

    fig.set_size_inches(16, 14)
    fig.tight_layout()
    plt.savefig(fname=output, dpi=dpi)


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'files', nargs='*',
        help='Radius of Gyration output from LAMMPS.')
    custom.add_argument(
        '--dpi', type=int, default=600,
        help='Resolution in dots per inch (default: %(default)s)')
    custom.add_argument(
        '--confidence', type=CI, default=0.95,
        help='Size of confidence interval to plot (default: %(default)s)')
    custom.add_argument(
        '--format', type=str,
        help='Output file format (default: auto)')
    requiredNamed = custom.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '--output', required=True,
        help='Output filename.')
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


def CI(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f'{x} not a floating-point literal')
    if x <= 0.0 or x >= 1.0:
        raise argparse.ArgumentTypeError(f'{x} not in range (0.0, 1.0)')
    return x


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
