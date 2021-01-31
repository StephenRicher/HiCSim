#!/usr/bin/env python3

""" Read lammps input dat file and define TAD status of each bead position  """


import sys
import json
import logging
import argparse
import fileinput
from argUtils import setDefaults, createMainParent


__version__ = '1.0.0'


def extractTADboundaries(dat: str):
    nBeads, CTCFboundaries =  parseInput(dat)
    outTADid = -1
    inTADid = 1
    beadStatus = {}
    for bead in range(1, nBeads + 1):
        status = tadStatus(bead, CTCFboundaries)
        if status == 'boundary':
            section = 0
            # Dont modify ID twice for consequtive boundaries
            if prevStatus != 'boundary':
                outTADid -= 1
                inTADid += 1
        elif status == 'inTAD':
            section = inTADid
        else:
            section = outTADid
        prevStatus = status
        beadStatus[bead] = section
    json.dump(beadStatus, sys.stdout)


def tadStatus(beadID, CTCFboundaries):
    """ Check if bead ID is inside a TAD, outside or a boundary """
    for forwardCTCF, reverseCTCF in CTCFboundaries:
        if (beadID == forwardCTCF) or (beadID == reverseCTCF):
            return 'boundary'
        elif (beadID > forwardCTCF) and (beadID < reverseCTCF):
            return 'inTAD'
    return 'outTAD'


def parseInput(dat: str):
    """ Extract nBeads and CTCF boundary positions """
    sections = ['Atoms', 'Bonds']
    nBeads = 0
    CTCFboundaries = []
    currentSection = None
    with fileinput.input(dat) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line in sections:
                currentSection = line
            elif currentSection == 'Atoms':
                if '# DNA ' in line:
                    nBeads += 1
            elif currentSection == 'Bonds':
                if line.endswith('# CTCF Bond'):
                    forwardCTCF, reverseCTCF = line.split()[2:4]
                    CTCFboundaries.append([int(forwardCTCF), int(reverseCTCF)])
    return nBeads, CTCFboundaries


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=extractTADboundaries)
    parser.add_argument(
        'dat', nargs='?', help='Lammps input dat file (default: stdin)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
