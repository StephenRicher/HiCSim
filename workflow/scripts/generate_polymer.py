#!/usr/bin/env python3

""" Script to generate a random linear polymer as an input file for LAMMPS """

import sys
import math
import random
import argparse
import numpy as np
import pyCommonTools as pct
from utilities import commaPair
from collections import namedtuple
from collections import defaultdict

#https://stackoverflow.com/questions/5154716/using-argparse-to-parse-arguments-of-form-arg-val
def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    base_args = pct.get_base_args()
    subparser = pct.make_subparser(parser)

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
        '--seed', default=42, type=int,
        help='Seed for random number generator.')

    random_subparser = subparser.add_parser(
        'random',
        description=random_linear_polymer.__doc__,
        help=random_linear_polymer.__doc__,
        parents=[base_args, box_sizes, seed_arg],
        epilog=parser.epilog)
    random_subparser.add_argument(
        '--n_molecules', default=1000, type=int,
        help='Number of molecules in polymer.')
    random_subparser.add_argument(
        '--n_types', default=4, type=int,
        help='Number of atom types in polymer.')
    random_subparser.add_argument(
        '--n_clusters', default=10, type=int,
        help='Number of clusters in poymers.')
    random_subparser.set_defaults(function=random_linear_polymer)

    track_subparser = subparser.add_parser(
        'model',
        description=create_polymer.__doc__,
        help=create_polymer.__doc__,
        parents=[base_args, box_sizes, seed_arg],
        epilog=parser.epilog)
    track_subparser.add_argument(
        'sequence',)
    track_subparser.add_argument(
        '--monomers', metavar='NBEADS,TYPE', type=commaPair,
        action='append', default=[],
        help='Add NBEADS number of monomers of TYPE.'
        'Call multiple times to add multiple monomer types.')
    track_subparser.add_argument(
        '--ctcf', default=False, action='store_true',
        help='Process beads labelled F and R as CTCF sites. Each bead is given '
        'a unique ID and convergent orientation pair coeffs are output.')
    track_subparser.add_argument(
        '--pairCoeffs',
        help='File of atom type coefficient pairings.')
    track_subparser.add_argument(
        '--coeffOut', default=None,
        help='Outfile to write atom pair coefficients (default: stderr)')
    track_subparser.add_argument(
        '--groupOut', default=None,
        help='Outfile to write atom group ID assignments (default: stderr)')
    track_subparser.add_argument(
        '--basesPerBead', required=True, type=int,
        help='Number of bases used to represent 1 bead.')
    track_subparser.set_defaults(function=create_polymer)

    # Assert correct version - required for preserving dictionary insert order
    assert sys.version_info >= (3, 7)
    return (pct.execute(parser))


class Sequence:

    def __init__(self, sequence, basesPerBead,
            name='DNA', ctcf=False, beadID=1):
        self._beadID = beadID
        self.basesPerBead = basesPerBead
        # Should beads set as 'F' or 'R' be processed as CTCF?
        self._ctcf = ctcf
        self.name = name
        self.load(sequence)

    def load(self, sequence):
        self.sequence = {}
        self.types = []
        with open(sequence) as f:
            for line in f:
                type = line.strip()
                if self._ctcf and type in ['F', 'R', 'B']:
                    type = f'{type}-{self._beadID}'
                if type not in self.types:
                    self.types.append(type)
                self.sequence[self._beadID] = type
                self._beadID += 1

    @property
    def nBeads(self):
        return len(self.sequence)

    @property
    def nAngles(self):
        return self.nBeads - 2

    @property
    def nBonds(self):
        return self.nBeads - 1

    @property
    def nTypes(self):
        return len(self.types)


class lammps:

    def __init__(self):
        self._log = pct.create_logger()
        self.sequences = []
        self.monomers = []
        self._beadID = 1
        self._typeID = 1
        self.typeIDs = {}

    def loadSequence(self, sequence_file, basesPerBead, ctcf):
        sequence = Sequence(sequence_file, basesPerBead=basesPerBead,
            ctcf=ctcf, beadID=self._beadID)
        self.sequences.append(sequence)
        self._beadID += sequence.nBeads
        for type in sequence.types:
            self.addType(type)

    def loadMonomer(self, type, nbeads):
        monomer = {}
        self.addType(type)
        for bead in range(nbeads):
            monomer[self._beadID] = type
            self._beadID += 1
        self.monomers.append(monomer)

    def loadBox(self, xlo, xhi, ylo, yhi, zlo, zhi):
        Dim = namedtuple('Dim', 'lo hi')
        Box = namedtuple('Box', ['x', 'y', 'z'])
        self.box = Box(Dim(xlo, xhi), Dim(ylo, yhi), Dim(zlo, zhi))

    def addType(self, type):
        if type not in self.typeIDs:
            self.typeIDs[type] = self._typeID
            self._typeID += 1

    @property
    def nBeads(self):
        n = sum([sequence.nBeads for sequence in self.sequences])
        n += sum([len(monomer) for monomer in self.monomers])
        return n

    @property
    def nAngles(self):
        return sum([sequence.nAngles for sequence in self.sequences])

    @property
    def nBonds(self):
        return sum([sequence.nBonds for sequence in self.sequences])

    @property
    def nTypes(self):
        return len(self.typeIDs)

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
            f'1 bond types\n'
            f'1 angle types\n\n'
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

    def writePolymers_old(self, r=1.1):
        for sequence in self.sequences:
            angle = 0
            for i, (beadID, type) in enumerate(sequence.sequence.items()):
                if i == 0:
                    x = random.random()
                    y = random.random()
                    z = random.random()
                else:
                    phi = random.random() * 2 * math.pi
                    theta = random.random() * math.pi
                    x = prev_x + (r * math.sin(theta) * math.cos(phi))
                    y = prev_y + (r * math.sin(theta) * math.sin(phi))
                    z = prev_z + (r * math.cos(theta))
                typeID = self.typeIDs[type]
                sys.stdout.write(f'{beadID} 1 {typeID} {x} {y} {z} 0 0 0\n')
                prev_x = x
                prev_y = y
                prev_z = z

    def writePolymers(self):
        for sequence in self.sequences:
            # Randomly select rosette length between 30kbp - 100kbp
            rosette_length = random.randint(30,100) * 1000
            beads_per_loop = rosette_length / sequence.basesPerBead
            # Randomly select number of loops per turn between 4 and 12
            loops_per_turn = random.randint(4,12)
            beads_per_turn = beads_per_loop * loops_per_turn
            n_turns = sequence.nBeads / beads_per_turn

            k = loops_per_turn / 2
            vx = 0.38
            theta = np.linspace(0, n_turns * 2 * np.pi, sequence.nBeads)
            x = 12 * (vx + (1 - vx) * (np.cos(k * theta))**2 * np.cos(theta))
            y = 12 * (vx + (1 - vx) * (np.cos(k * theta))**2 * np.sin(theta))
            z = theta / (2 * np.pi)

            for i, (beadID, type) in enumerate(sequence.sequence.items()):
                typeID = self.typeIDs[type]
                sys.stdout.write(f'{beadID} 1 {typeID} {x[i]} {y[i]} {z[i]} 0 0 0\n')

    def writeMonomers(self):
        for monomer in self.monomers:
            for beadID, type in monomer.items():
                x = random.uniform(self.box.x.lo, self.box.x.hi)
                y = random.uniform(self.box.y.lo, self.box.y.hi)
                z = random.uniform(self.box.z.lo, self.box.z.hi)
                typeID = self.typeIDs[type]
                sys.stdout.write(f'{beadID} 1 {typeID} {x} {y} {z} 0 0 0\n')


    def writeBonds(self):
        sys.stdout.write('Bonds\n\n')
        for sequence in self.sequences:
            nBonds = sequence.nBonds
            for i, beadID in enumerate(sequence.sequence):
                if i == sequence.nBonds:
                    break
                sys.stdout.write(f'{beadID} 1 {beadID} {beadID+1}\n')
        sys.stdout.write('\n')


    def writeAngles(self):
        sys.stdout.write('Angles\n\n')
        for sequence in self.sequences:
            nAngles = sequence.nAngles
            for i, beadID in enumerate(sequence.sequence):
                if i == sequence.nAngles:
                    break
                sys.stdout.write(f'{beadID} 1 {beadID} {beadID+1} {beadID+2}\n')
        sys.stdout.write('\n')


    def detectConvertCTCF(self, sequence):
        """ Detect valid CTCF convergent sites. For each interval, find the nearest
            F bead to the left and the nearest R bead to the right. This models the
            loop extrusion hypothesis.
        """
        pairs = set()
        typeList = sequence.types
        for beadIdx in range(len(typeList) - 1):
            # Reverse down the list until reach forward CTCF
            for forwardBead in reversed(typeList[:beadIdx + 1]):
                if forwardBead.startswith(('F', 'B')):
                    forwardTypeID = self.typeIDs[forwardBead]
                    break
            else:
                continue
            # Move up the list until reach reverse CTCF
            for reverseBead in typeList[beadIdx + 1:]:
                if reverseBead.startswith(('R', 'B')):
                    reverseTypeID = self.typeIDs[reverseBead]
                    break
            else:
                continue
            pairs.add((forwardTypeID, reverseTypeID))
        return pairs


    def writeConvergentCTCFs(self, coeff, fh=sys.stderr):
        for sequence in self.sequences:
            for forward, rev in self.detectConvertCTCF(sequence):
                fh.write(f'pair_coeff {forward} {rev} {coeff}\n')


    def writeCoeffs(self, type1, type2, coeff, fh=sys.stderr):
        try:
            type1ID = self.typeIDs[type1]
            type2ID = self.typeIDs[type2]
            # Ensure the smaller type_ID is written first
            if type1ID > type2ID:
                temp = type1ID
                type1ID = type2ID
                type2ID = temp
            fh.write(f'pair_coeff {type1ID} {type2ID} {coeff}\n')
        except KeyError:
            pass

    def writeGroup(self, fh=sys.stderr):
        fh.write(
            '####################################\n'
            '########       GROUPS        #######\n'
            '####################################\n')
        for sequence in self.sequences:
            firstBeadID = min(sequence.sequence)
            lastBeadID = max(sequence.sequence)
            name = sequence.name
            fh.write(f'group {name} id {firstBeadID}:{lastBeadID}\n')



def create_polymer(sequence, seed, monomers, ctcf, basesPerBead,
        pairCoeffs, coeffOut, groupOut, xlo, xhi, ylo, yhi, zlo, zhi):

    random.seed(seed)
    dat = lammps()
    dat.loadBox(xlo, xhi, ylo, yhi, zlo, zhi)
    dat.loadSequence(sequence, basesPerBead, ctcf)
    for monomer in monomers:
        nbeads = int(monomer[0])
        type = monomer[1]
        dat.loadMonomer(type, nbeads)
    dat.writeLammps()
    with pct.open(pairCoeffs) as f:
        with pct.open(coeffOut, 'w') as out:
            for line in f:
                line = line.strip().split()
                coeff = ' '.join(line[2:])
                atom1 = line[0]
                atom2 = line[1]
                # IF F-R or R-C (CTCF interaction)
                if sorted([atom1, atom2]) == ['F', 'R']:
                    dat.writeConvergentCTCFs(coeff, fh=out)
                else:
                    dat.writeCoeffs(atom1, atom2, coeff, fh=out)
    with pct.open(groupOut, 'w') as out:
        dat.writeGroup(fh=out)


def random_linear_polymer(
        n_molecules, n_types, n_clusters, seed, xlo, xhi, ylo, yhi, zlo, zhi):

    """ Generate a random linear polymer as an input file for LAMMPS """

    random.seed(seed)

    write_header(n_molecules, n_types, xlo, xhi, ylo, yhi, zlo, zhi)
    type_list = list(range(1, n_types + 1))
    write_masses(type_list)

    sys.stdout.write('\nAtoms\n\n')
    r = 1.1
    cluster_size = int(n_molecules / n_clusters)
    type = random.choice(type_list)
    type_boundary_idx = 0  # Set initial index of type boundary change
    for n in range(1, n_molecules+1):

        if n > type_boundary_idx + cluster_size:
            type_boundary_idx = n
            # Generate list of types excluding the previous type used
            type_choices = [n for n in type_list if n != type]
            type = random.choice(type_choices)

        if n == 1:
            x = random.random()
            y = random.random()
            z = random.random()
        else:
            phi = random.random() * 2 * math.pi
            theta = random.random() * math.pi
            x = prev_x + (r * math.sin(theta) * math.cos(phi))
            y = prev_y + (r * math.sin(theta) * math.sin(phi))
            z = prev_z + (r * math.cos(theta))

        sys.stdout.write(f'{n} 1 {type} {x} {y} {z} 0 0 0\n')
        prev_x = x
        prev_y = y
        prev_z = z
    sys.stdout.write('\n')
    write_bonds(n_molecules)

    write_angles(n_molecules)


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


if __name__ == '__main__':
    sys.exit(main())
