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
        infile: str, timestep: int,  simTime: int, dat: str,
        softWarmUp: int, harmonicWarmUp: int,
        restartPrefix: str, restartStep: int, ctcfSteps: int, seed: int,
        warmUpOut: str, simOut: str, groups: str, pairCoeff: str,
        TFswap: int, cosinePotential: float, radiusGyrationOut: str) -> None:

    if (not restartPrefix) and (restartStep != 0):
        logging.warning('Restart ignored as restartPrefix not provided.')
        restartStep = 0
    elif restartPrefix and restartStep == 0:
        logging.warning('Restart ignored as restartStep set to 0.')
        restartPrefix = ''

    if hasCTCFbond(dat):
        ctcfBond = True
        harmonicWarmUp = round(harmonicWarmUp / ctcfSteps)
        if harmonicWarmUp == 0:
            logging.error(
                'Progressive harmonic warm up cannot be less than 1. '
                'Reduce --ctcfSteps or increase --harmonicWarmup.')
            return 1
    else:
        # If no CTCF add any harmonicWarmUp to softWarmUp
        softWarmUp += harmonicWarmUp
        harmonicWarmUp = 0
    startTimestep = -(softWarmUp + (harmonicWarmUp * ctcfSteps) + timestep)


    with fileinput.input(infile) as fh:
        for line in fh:
            line = line.strip()
            if line == '## PAIR_COEFF ##':
                catFile(pairCoeff)
            elif line == '## GROUPS ##':
                catFile(groups)
            elif line.startswith('# CTCF BOND #') and ctcfBond:
                print(re.sub('# CTCF BOND # ', '', line))
            elif line.startswith('# CTCF Progressive #') and ctcfBond:
                harmonicCoeffs = np.linspace(0, 4, ctcfSteps)
                template = next(fh).strip('# \n')
                for coeff in harmonicCoeffs:
                    print(re.sub('\${harmonicCoeff}', str(coeff), template))
                    print(f'run {harmonicWarmUp}')
            else:
                line = re.sub('\${seed}', str(seed), line)
                line = re.sub('\${simTime}', str(simTime), line)
                line = re.sub('\${timestep}', str(timestep), line)
                line = re.sub('\${softWarmUp}', str(softWarmUp), line)
                line = re.sub('\${startTimestep}', str(startTimestep), line)
                line = re.sub('\${cosinePotential}', str(cosinePotential), line)
                line = re.sub('\${restartStep}', str(restartStep), line)
                line = re.sub('\${TFswap}', str(TFswap), line)
                line = re.sub('\${restartPrefix}', restartPrefix, line)
                line = re.sub('\${infile}', dat, line)
                line = re.sub('\${radiusGyrationOut}', radiusGyrationOut, line)
                line = re.sub('\${warmUpOut}', warmUpOut, line)
                line = re.sub('\${simOut}', simOut, line)
                print(line)


def hasCTCFbond(dat):
    """ Check if dat file contains a CTCF bond """
    with open(dat) as fh:
        for line in fh:
            if line.endswith('bond types\n'):
                nBonds = int(line.split()[0])
                return nBonds == 2


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
        '--timestep', type=int, default=500,
        help='Interval to write simulation state (default: %(default)s)')
    parser.add_argument(
        '--softWarmUp', type=int, default=100000,
        help='Warm up time with soft interactions (default: %(default)s)')
    parser.add_argument(
        '--harmonicWarmUp', type=int, default=100000,
        help='Warm up time for increasing harmonic bond (default: %(default)s)')
    parser.add_argument(
        '--simTime', type=int, default=100000,
        help='Total simulation time after warm up (default: %(default)s)')
    parser.add_argument(
        '--restartPrefix', default='',
        help='File prefix for restart files (default: no restart)')
    parser.add_argument(
        '--restartStep', type=int, default=0,
        help='Frequency to write restart state (default: %(default)s)')
    parser.add_argument(
        '--ctcfSteps', type=int, default=100,
        help='Steps in CTCF bond strengh increment (default: %(default)s)')
    parser.add_argument(
        '--TFswap', type=int, default=10000,
        help='Swap frequency of TF active/inactivation (default: %(default)s)')
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
