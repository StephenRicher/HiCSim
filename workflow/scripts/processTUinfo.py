#!/usr/bin/env python3

""" Determine TU activation state at each time point """

import sys
import logging
import argparse
import numpy as np
import pandas as pd
from argUtils import setDefaults, createMainParent
from utilities import getAtomCount, readJSON
from scipy.spatial.distance import cdist

__version__ = '1.0.0'


def processTUinfo(atomGroups: str, TADstatus: str, sim: str,
                  distance: float, out: str):

    # Read atomGroup and retrieve TU atom indexes
    atomGroups = readJSON(atomGroups)
    if ('TU' not in atomGroups) or ('TF' not in atomGroups):
        logging.error('No TU-TF pairs to process in simulation.')
        return 1

    atomCount = getAtomCount(atomGroups)

    # Read simulation file 1 timestep at a time
    info = pd.read_csv(
        sim, chunksize=atomCount,
        usecols=['timestep', 'ID', 'type', 'x', 'y', 'z'])

    allTUinfo = []
    for sim in info:
        # Extract data of only TUs
        TUdata = sim.loc[sim['ID'].isin(atomGroups['TU'])]
        TUs = TUdata[['x', 'y', 'z']].to_numpy()

        # Read TFs and exclude type 2 (inactive TF)
        TFas = sim.loc[(sim['ID'].isin(atomGroups['TF']))
                       & (sim['type'] != '2'), ['x', 'y', 'z']].to_numpy()

        # Compute distance of active TFs to each TU bead
        result = cdist(TUs, TFas, 'euclidean')

        TUdata = TUdata.copy()
        TUdata['TFnear'] = np.amin(result, axis=1)
        TUdata['active'] = TUdata['TFnear'] < distance
        allTUinfo.append(TUdata)

    allTUinfo = pd.concat(allTUinfo)
    # Read per-bead/timepoint TAD status
    TADstatus = pd.read_csv(TADstatus)
    # Exclude non-TU beads
    TADstatus = TADstatus[TADstatus['ID'].isin(atomGroups['TU'])]
    allTUinfo = pd.merge(allTUinfo, TADstatus)
    allTUinfo.to_csv(out, index=False)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument(
        'atomGroups', help='Atom group assignments in JSON format.')
    parser.add_argument(
        'TADstatus', help='Per bead/timepoint TAD status.')
    parser.add_argument(
        'sim', help='Input custom lammps simulation file.')
    parser.add_argument(
        '--distance', default=1.8, type=float,
        help='Distance threshold for transcription (default: %(default)s)')
    parser.add_argument(
        '--out', default=sys.stdout, help='File to TU info (default: stdout)')
    parser.set_defaults(function=processTUinfo)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
