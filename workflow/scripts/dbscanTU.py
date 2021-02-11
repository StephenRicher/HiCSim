#!/usr/bin/env python3

""" Cluster TUs with DBSCAN """


import sys
import logging
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from itertools import permutations
from sklearn.cluster import DBSCAN
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform


__version__ = '1.0.0'


def main(infile: str, out: str, outplot: str, eps: float, minSamples: float) -> None:

    fullSim = pd.read_csv(infile)
    TUids = fullSim['ID'].unique()
    fullSim = fullSim.groupby('timestep')
    clusterLabels = []
    for time, timestep in fullSim:
        positions = timestep[['x','y','z']].to_numpy()
        distances = squareform(pdist(positions, 'euclidean'))
        clusterLabels.extend(
            DBSCAN(eps=6, min_samples=3).fit(distances).labels_)

    # Ungroup and add labels column
    fullSim = fullSim.head(fullSim.ngroup().size)
    fullSim['labels'] = clusterLabels

    # Set inactive TUs to no cluster
    fullSim.loc[fullSim.active == False, 'labels'] = -1

    # Filter time periods with no defined clusters
    fullSim = fullSim.groupby('timestep').filter(lambda x: (x['labels'] != -1).any())

    if outplot:
        clusterLabels = fullSim[fullSim.labels != -1].labels.unique()
        nClusters = len(clusterLabels)
        cmap = sns.color_palette("Dark2", nClusters)

        clusterGroups = fullSim.pivot(
            index='ID', columns='timestep', values='labels')
        if clusterGroups.empty:
            logging.error('No valid clusters - creating empty plot file.')
            Path(outplot).touch()
        else:
            fig, ax = plt.subplots()
            ax = sns.heatmap(clusterGroups, ax=ax, mask=(clusterGroups==-1), cmap=cmap)
            colorbar = ax.collections[0].colorbar
            r = colorbar.vmax - colorbar.vmin
            colorbar.set_ticks(
                [colorbar.vmin + r / nClusters * (0.5 + i) for i in range(nClusters)])
            colorbar.set_ticklabels(clusterLabels)
            fig.savefig(outplot, dpi=300, bbox_inches='tight')

    clusterPairs = defaultdict(lambda: defaultdict(int))
    # Initialise all pairs as zero
    for id1, id2 in permutations(TUids, 2):
        clusterPairs[id1][id2] = 0

    # Count each time pairs occur in cluster
    for (time, label), group in fullSim.groupby(['timestep', 'labels']):
        if label == -1:
            continue
        ids = sorted(list(group.ID))
        for id1, id2 in permutations(group.ID, 2):
            clusterPairs[id1][id2] += 1

    clusterPairs = pd.concat(
        {k: pd.DataFrame(v, index=[0]).T for k, v in clusterPairs.items()}, axis=0)
    clusterPairs = clusterPairs.reset_index()
    clusterPairs.columns = ['TU1', 'TU2', 'count']
    timesteps = len(fullSim.timestep.unique())
    clusterPairs['proportion'] =  clusterPairs['count'] / timesteps
    clusterPairs.to_csv(out, index=False)


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument('infile', default=sys.stdin,
        help='Input TU stats.')
    parser.add_argument('--out', default=sys.stdout,
        help='File to save TU-TU cluster information.')
    parser.add_argument('--outplot',
        help='File to save TU cluster plot.')
    parser.add_argument(
        '--eps', type=float, default=6,
        help='EPS parameter for DBSCAN (default: %(default)s)')
    parser.add_argument(
        '--minSamples', type=float, default=2,
        help='min_samples parameter for DBSCAN (default: %(default)s)')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
