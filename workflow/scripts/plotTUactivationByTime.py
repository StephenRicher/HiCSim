#!/usr/bin/env python3

""" Compute TU-TU distance per time point and plot """


import sys
import logging
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import product
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform


__version__ = '1.0.0'


def main(file: str, outdir: str, **kwargs) -> None:

    fullSim = pd.read_csv(file)

    TUids = fullSim['id'].unique()
    fullSim = fullSim.groupby('time')
    TU1, TU2 = generatePairwise(TUids)

    allTUs = []
    for time, timestep in fullSim:
        refActive, compActive = generatePairwise(timestep['active'])
        positions = timestep[['x','y','z']].to_numpy()
        compActive = [1 if x else -1 for x in compActive]
        # Convert to squareform (to get self comparison)
        distances = squareform(pdist(positions, 'euclidean'))
        distances = distances.reshape(len(TU1))
        allTUs.append(pd.DataFrame(
            {'refTU'    : TU1                ,
             'compTU'   : TU2                ,
             'time'     : [time] * len(TU1)  ,
             'dist'     : distances * compActive   ,}))
    allTUs = pd.concat(allTUs)

    for TU in TUids:
        refTU = allTUs[allTUs['refTU'] == TU].pivot(
            index='compTU', columns='time', values='dist')

        # Compute larger number reflect closer association
        refTU = (1/refTU)
        # Get maximum and minimum non-infinite numbers
        maxNonInf = refTU.replace([np.inf, -np.inf], np.nan).max().max()
        minNonInf = refTU.replace([np.inf, -np.inf], np.nan).min().min()

        refTU = refTU.replace([np.inf, -np.inf], [maxNonInf, minNonInf])

        fig, ax = plt.subplots()
        ax = sns.heatmap(refTU, cmap='bwr', center=0, ax=ax)
        fig.tight_layout()
        fig.savefig(f'{outdir}/{TU}-activity.png', dpi=300, bbox_inches='tight')

def generatePairwise(names):
    names1 = []
    names2 = []
    names = product(names, repeat=2)
    for name1, name2 in names:
        names1.append(name1)
        names2.append(name2)
    return names1, names2


def parse_arguments():

    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument('file',
        help='Input TU activation table.')
    custom.add_argument('--outdir', default='.',
        help='Directory to save output plots.')


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
