#!/usr/bin/env python3

""" Script to generate a random linear polymer as an input file for LAMMPS """

import sys
import math
import random
import argparse
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

    beads = argparse.ArgumentParser(add_help=False)
    beads.add_argument(
        '--single_beads', metavar='NBEADS,TYPE', type=commaPair,
        action='append', default=[],
        help='Add NBEADS number of isolated beads of TYPE.'
        'Call multiple times to add multiple beads.')

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
        parents=[base_args, box_sizes, seed_arg, beads],
        epilog=parser.epilog)
    track_subparser.add_argument(
        'sequence',)
    track_subparser.add_argument(
        '--bead_pair', nargs='+',
        metavar="KEY=VALUE", action=ParseDict,
        help='Set bead pair interactions. e.g. For 2 beads of type \'I\''
        'and \'J\' we may set the interaction using I-J=1.0,1.0,2.5' )
    track_subparser.add_argument(
        '--ctcf', default=False, action='store_true',
        help='Process beads labelled F and R as CTCF sites. Each bead is given '
        'a unique ID and convergent orientation pair coeffs are output.')
    track_subparser.add_argument(
        '--ctcf_coeff', default='1.5,1.0,1.8',
        help='Define pairing coefficient between convergent CTCF sites.')
    track_subparser.add_argument(
        '--ctcf_out',
        help='Outfile to write convergent ctcf pair_coeffs (default: stderr)')
    track_subparser.set_defaults(function=create_polymer)

    # Assert correct version - required for preserving dictionary insert order
    assert sys.version_info >= (3, 7)
    return (pct.execute(parser))



class polymer:

    def __init__(self, sequence, single_beads, ctcf,
            xlo, xhi, ylo, yhi, zlo, zhi):
        log = pct.create_logger()
        Dim = namedtuple('Dim', 'lo hi')
        Box = namedtuple('Box', ['x', 'y', 'z'])
        self.box = Box(Dim(xlo, xhi), Dim(ylo, yhi), Dim(zlo, zhi))
        self.ctcf = []
        self.boundAtoms = defaultdict(int)
        self.sequence = []
        with open(sequence) as f:
            for line in f:
                bead = line.strip()
                self.sequence.append(bead)
                if ctcf and bead in ['F', 'R']:
                    id = len(self.ctcf) + 1
                    bead = f'{bead}-{id}'
                    self.ctcf.append(bead)
                self.boundAtoms[bead] += 1

        self.unboundAtoms = defaultdict(int)
        for single_bead in single_beads:
            nbeads = int(single_bead[0])
            bead = single_bead[1]
            if bead in self.boundAtoms:
                log.error(
                    f'Single bead type {bead} is alread used in {sequence}')
                sys.exit(1)
            elif bead in ['F', 'R']:
                log.error(
                    'Single bead type cannot be F or R.')
                sys.exit(1)
            else:
                self.unboundAtoms[bead] = int(nbeads)

    @property
    def nTotalAtoms(self):
        return self.nBoundAtoms + self.nUnboundAtoms

    @property
    def nBoundAtoms(self):
        return sum(self.boundAtoms.values())

    @property
    def nUnboundAtoms(self):
        return sum(self.unboundAtoms.values())

    @property
    def nAngles(self):
        return self.nBoundAtoms - 2

    @property
    def nBonds(self):
        return self.nBoundAtoms - 1

    @property
    def nTypes(self):
        return len(self.boundAtoms) + len(self.unboundAtoms)

    @property
    def typeList(self):
        return list(self.boundAtoms) + list(self.unboundAtoms)

    def typeNumber(self, type):
        return self.typeList.index(type) + 1

    def writeHeader(self):
        sys.stdout.write(
            'LAMMPS data file from restart file: timestep = 0, procs = 1\n\n'
            f'{self.nTotalAtoms} atoms\n'
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
        for n, type in enumerate(self.typeList):
            sys.stdout.write(f'{n+1} 1 # {type}\n')
        sys.stdout.write('\n')

    def writeAtoms(self):
        sys.stdout.write('Atoms\n\n')
        self.writeBoundAtoms()
        self.writeUnboundAtoms()
        sys.stdout.write('\n')

    def writeBoundAtoms(self, r=1.1):
        ctcf_beads = self.ctcf
        for i, type in enumerate(self.sequence):
            if self.ctcf and type in ['F', 'R']:
                type = ctcf_beads[0]
                ctcf_beads = ctcf_beads[1:]
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
            type_number = self.typeNumber(type)
            sys.stdout.write(f'{i+1} 1 {type_number} {x} {y} {z} 0 0 0\n')
            prev_x = x
            prev_y = y
            prev_z = z

    def writeUnboundAtoms(self):
        atom_id = self.nBoundAtoms + 1
        for type, nbeads in self.unboundAtoms.items():
            type_number = self.typeNumber(type)
            for n in range(nbeads):
                x = random.random()
                y = random.random()
                z = random.random()
                sys.stdout.write(
                    f'{atom_id} 1 {type_number} {x} {y} {z} 0 0 0\n')
                atom_id += 1

    def writeBonds(self):
        sys.stdout.write('Bonds\n\n')
        for b in range(1, self.nBonds + 1):
            sys.stdout.write(f'{b} 1 {b} {(b%self.nBoundAtoms)+1}\n')
        sys.stdout.write('\n')


    def writeAngles(self):
        sys.stdout.write('Angles\n\n')
        for a in range(1, self.nAngles + 1):
            sys.stdout.write(f'{a} 1 {a} {a+1} {a+2}\n')
        sys.stdout.write('\n')



def create_polymer(sequence, seed, bead_pair, single_beads, ctcf,
        ctcf_coeff, ctcf_out, xlo, xhi, ylo, yhi, zlo, zhi):

    random.seed(seed)
    a = polymer(sequence, single_beads, ctcf, xlo, xhi, ylo, yhi, zlo, zhi)
    a.writeHeader()
    a.writeMasses()
    a.writeAtoms()
    a.writeBonds()
    a.writeAngles()
    exit(1)
    write_ctcf_interactions(type_list, ctcf_coeff, ctcf_out)

def detect_pairs(beads):
    """ Detect valid CTCF convergent sites. For each interval, find the nearest
        F bead to the left and the nearest R bead to the right. This models the
        loop extrusion hypothesis.
    """

    pairs = set()
    for bead_idx in range(0, len(beads) - 1):
        for forward_bead in reversed(beads[:bead_idx + 1]):
            if forward_bead.startswith('F'):
                forward_number = beads.index(forward_bead) + 1
                break
        else:
            continue
        for reverse_bead in beads[bead_idx + 1:]:
            if reverse_bead.startswith('R'):
                reverse_number = beads.index(reverse_bead) + 1
                break
        else:
            continue
        pairs.add((forward_number, reverse_number))
    return pairs


def write_ctcf_interactions(type_list, ctcf_coeff, ctcf_out):
    ctcf_coeff = ctcf_coeff.replace(',',' ')
    with pct.open(ctcf_out, mode='w', stderr=False) as f:
        for forward, reverse in detect_pairs(type_list):
            f.write(f'pair_coeff {forward} {reverse} {ctcf_coeff}\n')


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
