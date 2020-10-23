#!/usr/bin/env python3

""" Compute TU-TU distance per time point and plot """

import sys
import logging
import argparse
import sys
import json
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import seaborn as sns
import matplotlib.pyplot as plt

__version__ = '1.0.0'


def main(file: str, outdir: str, **kwargs) -> None:

    fullSim = pd.read_csv(file)

    TUids = fullSim['id'].unique()
    fullSim = fullSim.groupby('time')

    allTUs = []
    for time, timestep in fullSim:
        for TU in TUids:
            active = [1 if x else -1 for x in timestep['active']]
            referencePos = timestep.loc[timestep['id'] == TU, ['x','y','z']].to_numpy()
            comparePos = timestep[['x','y','z']].to_numpy()
            distances = cdist(referencePos, comparePos, 'euclidean')[0]
            timestep['cdist'] = distances * active
            allTUs.append(pd.DataFrame(
                {'refTU'    : [TU] * len(TUids)    ,
                 'compTU'   : TUids                ,
                 'time'     : [time] * len(TUids)  ,
                 'dist'     : distances * active   ,}))
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
