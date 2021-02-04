#!/usr/bin/env python3

""" Script to generate a random linear polymer as an input file for LAMMPS """

import sys
import math
import random
import logging
import argparse
import fileinput
import numpy as np
from collections import namedtuple
from utilities import readJSON
from collections import defaultdict

__version__ = '1.0.0'


class lammps:

    def __init__(self, nPolymerBeads, basesPerBead, beadTypes={}, nMonomers=0, polymerSeed=random.randint(1, 1e100),
            extrusion=False, monomerSeed=random.randint(1, 1e100),
            randomWalk=False):
        self.nPolymerBeads = nPolymerBeads
        self.basesPerBead = basesPerBead
        self._beadID = 1
        # Merge typeID mappings
        self.typeIDs = {**{'TFa': 1, 'TFi': 2, 'N': 3}, **beadTypes}
        self.randomWalk = randomWalk
        self._monomerSeed = monomerSeed
        self._polymerSeed = polymerSeed
        self.nMonomers = nMonomers
        self.extrusion = extrusion

    @property
    def nBeads(self):
        return self.nPolymerBeads + self.nMonomers

    @property
    def nAngles(self):
        return self.nPolymerBeads - 2

    @property
    def nBonds(self):
        return self.nPolymerBeads - 1

    @property
    def nTypes(self):
        return len(self.typeIDs)

    @property
    def nBondTypes(self):
        return 2 if self.extrusion else 1


    def loadBox(self, xlo, xhi, ylo, yhi, zlo, zhi):
        Dim = namedtuple('Dim', 'lo hi')
        Box = namedtuple('Box', ['x', 'y', 'z'])
        self.box = Box(Dim(xlo, xhi), Dim(ylo, yhi), Dim(zlo, zhi))

    def writeLammps(self):
        self.writeHeader()
        self.writeMasses()
        self.writeBeads()
        self.writeBonds()
        self.writeAngles()

    def writeHeader(self):
        sys.stdout.write(
            'LAMMPS data file from restart file: timestep = 0, procs = 1\n\n'
            f'{self.nBeads} atoms\n'
            f'{self.nBonds} bonds\n'
            f'{self.nAngles} angles\n\n'
            f'{self.nTypes} atom types\n'
            f'{self.nBondTypes} bond types\n'
            f'1 angle types\n'
            f'{"2 extra bond per atom" if self.extrusion else ""}\n'
            f'{self.box.x.lo} {self.box.x.hi} xlo xhi\n'
            f'{self.box.y.lo} {self.box.y.hi} ylo yhi\n'
            f'{self.box.z.lo} {self.box.z.hi} zlo zhi\n\n')


    def writeMasses(self):
        sys.stdout.write('Masses\n\n')
        for type, typeID in self.typeIDs.items():
            sys.stdout.write(f'{typeID} 1 # {type}\n')
        sys.stdout.write('\n')


    def writeBeads(self):
        sys.stdout.write('Atoms\n\n')
        self.writePolymers()
        self.writeMonomers()
        sys.stdout.write('\n')


    def makeRandomWalk(self, length, r=1.1):
        x = []
        y = []
        z = []
        for i in range(length):
            if i == 0:
                x.append(random.random())
                y.append(random.random())
                z.append(random.random())
            else:
                phi = random.random() * 2 * math.pi
                theta = random.random() * math.pi
                x.append(x[i-1] + (r * math.sin(theta) * math.cos(phi)))
                y.append(y[i-1] + (r * math.sin(theta) * math.sin(phi)))
                z.append(z[i-1] + (r * math.cos(theta)))
        return x, y, z


    def writePolymers(self):
        random.seed(self._polymerSeed)
        if self.randomWalk:
            x, y, z = self.makeRandomWalk(self.nPolymerBeads)
        else:
            # Randomly select rosette length between 30kbp - 100kbp
            rosette_length = random.randint(30,100) * 1000
            beads_per_loop = rosette_length / self.basesPerBead
            # Randomly select number of loops per turn between 4 and 12
            loops_per_turn = random.randint(4,12)
            beads_per_turn = beads_per_loop * loops_per_turn
            n_turns = self.nPolymerBeads / beads_per_turn

            k = loops_per_turn / 2
            vx = 0.38
            theta = np.linspace(0, n_turns * 2 * np.pi, self.nPolymerBeads)
            x = 12 * (vx + (1 - vx) * (np.cos(k * theta))**2 * np.cos(theta))
            y = 12 * (vx + (1 - vx) * (np.cos(k * theta))**2 * np.sin(theta))
            z = theta / (2 * np.pi)
        for i in range(self.nPolymerBeads):
            sys.stdout.write(f'{self._beadID} 1 3 {x[i]} {y[i]} {z[i]} 0 0 0 # DNA\n')
            self._beadID += 1


    def writeMonomers(self, prop=0.5):
        random.seed(self._monomerSeed)
        activeMonomers = int(self.nMonomers * prop)
        inactiveMonomers = self.nMonomers - activeMonomers
        monomers = (   activeMonomers   * ['TFa']
                     + inactiveMonomers * ['TFi'])
        for monomer in monomers:
            x = random.uniform(self.box.x.lo, self.box.x.hi)
            y = random.uniform(self.box.y.lo, self.box.y.hi)
            z = random.uniform(self.box.z.lo, self.box.z.hi)
            typeID = self.typeIDs[monomer]
            sys.stdout.write(f'{self._beadID} 1 {typeID} {x} {y} {z} 0 0 0 # TF\n')
            self._beadID += 1


    def writeBonds(self):
        sys.stdout.write('Bonds\n\n')
        bondID = 1
        for bead in range(self.nPolymerBeads):
            if bead == self.nBonds:
                break
            beadID = bead + 1
            sys.stdout.write(f'{bondID} 1 {beadID} {beadID + 1}\n')
            bondID += 1
        sys.stdout.write('\n')


    def writeAngles(self):
        sys.stdout.write('Angles\n\n')
        for bead in range(self.nPolymerBeads):
            if bead == self.nAngles:
                break
            beadID = bead + 1
            sys.stdout.write(f'{beadID} 1 {beadID} {beadID + 1} {beadID + 2}\n')
        sys.stdout.write('\n')


def main(nPolymerBeads: int, beadTypes: str, monomerSeed, polymerSeed, nMonomers, extrusion, basesPerBead, randomWalk,
        xlo, xhi, ylo, yhi, zlo, zhi, **kwargs):

    dat = lammps(
        nPolymerBeads=nPolymerBeads, basesPerBead=basesPerBead, beadTypes=readJSON(beadTypes),
        nMonomers=nMonomers, monomerSeed=monomerSeed, polymerSeed=polymerSeed,
        randomWalk=randomWalk, extrusion=extrusion)
    dat.loadBox(xlo, xhi, ylo, yhi, zlo, zhi)
    dat.writeLammps()


class ParseDict(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        d = {}
        if values:
            for item in values:
                split_items = item.split("=", 1)
                key = split_items[
                    0
                ].strip()  # we remove blanks around keys, as is logical
                value = split_items[1]

                d[key] = value

        setattr(namespace, self.dest, d)


def parse_arguments():

    box_sizes = argparse.ArgumentParser(add_help=False)
    box_sizes.add_argument(
        '--xlo', default=-50, type=float,
        help='Lower x-axis of simuluation box.')
    box_sizes.add_argument(
        '--xhi', default=50, type=float,
        help='Upper x-axis of simuluation box.')
    box_sizes.add_argument(
        '--ylo', default=-50, type=float,
        help='Lower y-axis of simuluation box.')
    box_sizes.add_argument(
        '--yhi', default=50, type=float,
        help='Upper y-axis of simuluation box.')
    box_sizes.add_argument(
        '--zlo', default=-50, type=float,
        help='Lower z-axis of simuluation box.')
    box_sizes.add_argument(
        '--zhi', default=50, type=float,
        help='Upper z-axis of simuluation box.')

    seed_arg = argparse.ArgumentParser(add_help=False)
    seed_arg.add_argument(
        '--polymerSeed', default=None, type=float,
        help='Seed for polymer structure.')
    seed_arg.add_argument(
        '--monomerSeed', default=None, type=float,
        help='Seed for monomer positions.')


    custom = argparse.ArgumentParser(add_help=False)
    custom.set_defaults(function=main)
    custom.add_argument(
        'nPolymerBeads', type=int, help='Number of atoms in polymer.')
    custom.add_argument(
        '--beadTypes', default=None,
        help='Additional bead type to type ID mappings in JSON format.')
    custom.add_argument(
        '--nMonomers', default=0, type=int,
        help='Set number of monomers for simulation.')
    custom.add_argument(
        '--extrusion', default=False, action='store_true',
        help='Process beads labelled F and R as CTCF sites. Each '
             'bead is given a unique ID and convergent orientation '
             'pair coeffs are output (default: %(default)s)')
    custom.add_argument(
        '--basesPerBead', required=True, type=int,
        help='Number of bases used to represent 1 bead.')
    custom.add_argument(
        '--randomWalk', default=False, action='store_true',
        help='Initialise polymer conformation using a random walk rather than '
             'a helicoidal conoformation  (default: %(default)s)')
    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        '--version', action='version', version=f'%(prog)s {__version__}')
    base.add_argument(
        '--verbose', action='store_const', const=logging.DEBUG,
        default=logging.INFO, help='verbose logging for debugging')

    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[base, custom, box_sizes, seed_arg])
    args = parser.parse_args()

    log_format='%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    logging.basicConfig(level=args.verbose, format=log_format)

    # Assert correct version - required for preserving dictionary insert order
    assert sys.version_info >= (3, 7)

    return args


if __name__ == '__main__':
    args = parse_arguments()
    return_code = args.function(**vars(args))
    logging.shutdown()
    sys.exit(return_code)
