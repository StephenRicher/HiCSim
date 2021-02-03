#!/usr/bin/env python3

""" Write explicit Lammps script """

import re
import sys
import random
import logging
import argparse
import fileinput
import numpy as np
from utilities import setDefaults

__version__ = '1.0.0'


def writeLammps(
        infile: str, timestep: float, writeInterval: float, dat: str,
        warmUp: float, restartPrefix: str, restartStep: int, extrusion: bool,
        seed: int, warmUpOut: str, simOut: str, groups: str,
        harmonicCoeff: float, pairCoeff: str, TFswap: float,
        cosinePotential: float, radiusGyrationOut: str) -> None:

    if (not restartPrefix) and (restartStep != 0):
        logging.warning('Restart ignored as restartPrefix not provided.')
        restartStep = 0
    elif restartPrefix and restartStep == 0:
        logging.warning('Restart ignored as restartStep set to 0.')
        restartPrefix = ''

    # Convert time unit to timesteps
    writeIntervalSteps = int(writeInterval / timestep)
    warmUpSteps = int(warmUp / timestep)
    TFswapSteps = int(TFswap / timestep)
    startTimestep = -(warmUpSteps + writeIntervalSteps)
    with fileinput.input(infile) as fh:
        for line in fh:
            line = line.strip()
            if line == '## PAIR_COEFF ##':
                catFile(pairCoeff)
            elif line == '## GROUPS ##':
                catFile(groups)
            else:
                if extrusion:
                    line = re.sub('# EXTRUSION # ', '', line)
                    line = re.sub('\${harmonicCoeff}', str(harmonicCoeff), line)
                line = re.sub('\${seed}', str(seed), line)
                line = re.sub('\${dt}', str(timestep), line)
                line = re.sub('\${writeInterval}', str(writeIntervalSteps), line)
                line = re.sub('\${warmUp}', str(warmUpSteps), line)
                line = re.sub('\${startTimestep}', str(startTimestep), line)
                line = re.sub('\${cosinePotential}', str(cosinePotential), line)
                line = re.sub('\${restartStep}', str(restartStep), line)
                line = re.sub('\${TFswap}', str(TFswapSteps), line)
                line = re.sub('\${restartPrefix}', restartPrefix, line)
                line = re.sub('\${infile}', dat, line)
                line = re.sub('\${radiusGyrationOut}', radiusGyrationOut, line)
                line = re.sub('\${warmUpOut}', warmUpOut, line)
                line = re.sub('\${simOut}', simOut, line)
                print(line)


def catFile(file):
    """ Write file contents to stdout """
    with open(file) as fh:
        for line in fh:
            print(line, end='')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'infile', nargs='?', default=[],
        help='Input Lammps template (default: stdin)')
    parser.add_argument(
        '--timestep', type=float, default=0.01,
        help='Size of timestep in time units (default: %(default)s)')
    parser.add_argument(
        '--writeInterval', type=float, default=10,
        help='Simulation status output interval in time units '
             '(default: %(default)s)')
    parser.add_argument(
        '--TFswap', type=float, default=100,
        help='Swap frequency of TF active/inactivation in time units '
             '(default: %(default)s)')
    parser.add_argument(
        '--warmUp', type=int, default=10000,
        help='Warm up time with soft interactions (default: %(default)s)')
    parser.add_argument(
        '--restartPrefix', default='',
        help='File prefix for restart files (default: no restart)')
    parser.add_argument(
        '--restartStep', type=int, default=0,
        help='Frequency to write restart state (default: %(default)s)')
    parser.add_argument(
        '--extrusion', default=False, action='store_true',
        help='Set to prepare script for loop extrusion (default: %(default)s)')
    parser.add_argument(
        '--harmonicCoeff', type=float, default=40.0,
        help='Harmonic bond strength for extrusion factors '
             '(default: %(default)s)')
    parser.add_argument(
        '--seed', type=int, default=random.randint(0, 10e9),
        help='Seed for simulation (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--dat', required=True, help='Path to Lammps dat file.')
    requiredNamed.add_argument(
        '--warmUpOut', required=True, help='Output file for warm up.')
    requiredNamed.add_argument(
        '--simOut', required=True, help='Output file for simulation.')
    requiredNamed.add_argument(
        '--groups', required=True, help='Atom ID group definitions.')
    requiredNamed.add_argument(
        '--pairCoeff', required=True,
        help='Atom type pair coefficient definitions.')
    requiredNamed.add_argument(
        '--cosinePotential', type=float, required=True,
        help='Set potential for cosine angle')
    requiredNamed.add_argument(
        '--radiusGyrationOut', required=True,
        help='Output file for radius of gyration data.')

    return setDefaults(parser, verbose=True, version=__version__)


if __name__ == '__main__':
    args = parseArgs()
    sys.exit(writeLammps(**vars(args)))
