#!/usr/bin/env python3

""" Plot DBSCAN clusters - mean proportion of time TU pairs are in cluster. """

import sys
import argparse
import numpy as np
import pandas as pd
from typing import List
import matplotlib.pyplot as plt
import seaborn as sns

def main(infiles: List, out: str) -> None:

    mergedDF = pd.concat(pd.read_csv(file) for file in infiles)

    mergedDF = mergedDF.groupby(['TU1', 'TU2']).mean().reset_index()

    # Set diagonal to NaN to hide trivial auto-correlation
    mergedDF.loc[mergedDF.TU1 == mergedDF.TU2, 'proportion'] = np.nan
    mergedDF = mergedDF.pivot(index='TU1', columns='TU2', values='proportion')

    # Flip vertically to ensure diagonal goes from bottom left to top right
    mergedDF = mergedDF.iloc[::-1]

    fig, ax = plt.subplots()
    ax = sns.heatmap(mergedDF, square=True, cmap='binary', vmin=0, ax=ax)
    # Ensure masked cells are not within 'bwr' colour map.
    ax.set_facecolor('xkcd:light grey')

    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches='tight')


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('out', help='Output heatmap.')
    parser.add_argument('infiles', nargs='*', help='Input DBSCAN TU-TU pairs.')


    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
