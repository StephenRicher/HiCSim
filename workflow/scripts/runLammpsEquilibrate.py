#!/usr/bin/env python3

import sys
import random
import logging
import argparse
from lammps import lammps
from argUtils import setDefaults, createMainParent


__version__ = '1.0.0'


def runEquilibrationLammps(
        inputDat: str, equilOut: str, timestep: float, writeInterval: float,
        equilTime: float, seed: int, equilInfo: str,
        cosinePotential: float, radiusGyrationOut: str) -> None:

    # Convert time unit to timesteps
    writeIntervalSteps = int(writeInterval / timestep)
    equilSteps = int(equilTime / timestep)
    startTimestep = -(equilSteps + writeIntervalSteps)

    global lmp
    lmp = lammps()
    lmp.command('units lj')
    lmp.command('atom_style angle')
    lmp.command('boundary p p p')
    lmp.command('neighbor 1.4 bin')
    lmp.command('neigh_modify every 1 delay 1 check yes')

    lmp.command(f'restart 0')

    lmp.command(f'read_data {inputDat} extra/special/per/atom 2')

    lmp.command('pair_style soft 1.12246')
    lmp.command('pair_modify shift yes')
    lmp.command('pair_coeff * * 20 1.12246')
    lmp.command('variable prefactor equal ramp(0,50)')
    lmp.command('fix pushapart all adapt 1 pair soft a * * v_prefactor')

    lmp.command('bond_style harmonic')
    lmp.command('bond_coeff 1 100.0 1.1')
    lmp.command('bond_coeff 2 0.0 1.1')
    lmp.command('bond_coeff 3 0.0 1.1')

    lmp.command('angle_style cosine')
    lmp.command(f'angle_coeff 1 {cosinePotential}')

    lmp.command('fix 1 all nve')
    lmp.command(f'fix 2 all langevin 1.0 1.0 1.0 {seed}')

    lmp.command('thermo 100')
    lmp.command('thermo_style custom step temp epair vol cpu')
    lmp.command(f'timestep {timestep}')
    lmp.command('reset_timestep 0')

    writeGroups(inputDat)

    lmp.command('compute RG DNA gyration')
    lmp.command('variable RG equal c_RG')

    lmp.command(f'variable time equal {startTimestep}+step')
    lmp.command(f'fix RG DNA print {writeIntervalSteps} "${{time}} ${{RG}}" file {radiusGyrationOut}')

    lmp.command(f'dump 1 all custom {writeIntervalSteps} {equilInfo} id type x y z ix iy iz')
    lmp.command("dump_modify 1 format line '%d %d %.5f %.5f %.5f %d %d %d' sort 1")

    lmp.command(f'run {equilSteps}')
    lmp.command('unfix pushapart')
    lmp.command('undump 1')
    lmp.command('unfix RG')

    # Write restart after equilibration
    lmp.command(f'write_restart {equilOut}')


def writeGroups(inputDat):
    """ Read lammps input and define bead groups which are polymer and TFs """
    nPolymerBeads = 0
    nMonomers = 0
    with open(inputDat) as fh:
        for line in fh:
            line = line.strip()
            if line.endswith('# DNA'):
                nPolymerBeads += 1
            elif line.endswith('# TF'):
                nMonomers += 1
    lmp.command(f'group DNA id 1:{nPolymerBeads}')
    lmp.command(f'group TF id {nPolymerBeads + 1}:{nPolymerBeads + 1 + nMonomers}')


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=runEquilibrationLammps)
    parser.add_argument('inputDat', help='Lammps input data.')
    parser.add_argument(
        'equilOut', help='Path to write equilibration Restart file')
    parser.add_argument(
        '--timestep', type=float, default=0.01,
        help='Size of timestep in time units (default: %(default)s)')
    parser.add_argument(
        '--writeInterval', type=float, default=10,
        help='Simulation status output interval in time units '
             '(default: %(default)s)')
    parser.add_argument(
        '--equilTime', type=int, default=10000,
        help='Equilibration time with soft interactions (default: %(default)s)')
    parser.add_argument(
        '--seed', type=int, default=random.randint(0, 10e9),
        help='Seed for simulation (default: %(default)s)')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--equilInfo', required=True, help='Output file for warm up.')
    requiredNamed.add_argument(
        '--cosinePotential', type=float, required=True,
        help='Set potential for cosine angle')
    requiredNamed.add_argument(
        '--radiusGyrationOut', required=True,
        help='Output file for radius of gyration data.')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
