#!/usr/bin/env python3

""" Compute TAD status for each TU-TU pair """


import sys
import argparse
import pandas as pd
from itertools import combinations
from argUtils import setDefaults, createMainParent

__version__ = '1.0.0'


def processTADstatus(TUinfo: str, out: str) -> None:

    TUinfo = pd.read_csv(TUinfo)

    # Initliase df of pairwise TUs for each timestep
    TUpairs = (TUinfo.groupby('timestep')['ID']
        .apply(lambda x: pd.DataFrame((i for i in combinations(x.values,2))))
        .reset_index().drop('level_1',axis=1)
        .rename(columns={0: 'TU1', 1: 'TU2'}))

    # Merge information for both TU1 and TU2
    TUpairs = pd.merge(
        TUpairs, TUinfo[['ID', 'timestep', 'TADstatus', 'TADgroup', 'active']],
        left_on=['timestep', 'TU1'],
        right_on=['timestep', 'ID']).drop(['ID'], axis=1)
    TUpairs.rename(
        columns={'TADstatus': 'TADstatus-TU1',
                 'TADgroup' : 'TADgroup-TU1',
                 'active'   : 'active-TU1'},
        inplace=True)
    TUpairs = pd.merge(
        TUpairs, TUinfo[['ID', 'timestep', 'TADstatus', 'TADgroup', 'active']],
        left_on=['timestep', 'TU2'],
        right_on=['timestep', 'ID']).drop(['ID'], axis=1)
    TUpairs.rename(
        columns={'TADstatus': 'TADstatus-TU2',
                 'TADgroup' : 'TADgroup-TU2',
                 'active'   : 'active-TU2'},
        inplace=True)
    TUpairs['TADstatus'] = TUpairs.apply(setPairStatus, axis=1)

    TUpairs.to_csv(out, index=False)


def setPairStatus(row):
    """ Label TU-TU pair status by relative TAD positions """
    seperation = abs(abs(row['TADgroup-TU1']) - abs(row['TADgroup-TU2']))
    TU1_inTAD = row['TADstatus-TU1'] > 0
    TU2_inTAD = row['TADstatus-TU2'] > 0
    if seperation == 0:
        if TU1_inTAD:
            return 'sameTAD'
        else:
            return 'sameNotTAD'
    elif seperation == 1:
        if TU1_inTAD == TU2_inTAD:
            if TU1_inTAD:
                return 'adjacentTAD'
            else:
                return 'adjacentNotTAD'
        else:
            return 'adjacentBoundary'
    else:
        return 'other'


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.add_argument('TUinfo', default=sys.stdin,
        help='Input TU activation table.')
    parser.add_argument('--out', default=sys.stdout,
        help='File to save TU TAD status metrics.')
    parser.set_defaults(function=processTADstatus)

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
