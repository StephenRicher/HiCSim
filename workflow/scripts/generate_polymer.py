#!/usr/bin/env python3

""" Script to generate a random linear polymer as an input file for LAMMPS """

import sys
import math
import random
import argparse
import pyCommonTools as pct
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


    return (pct.execute(parser))


def create_polymer(sequence, seed, bead_pair, ctcf, ctcf_coeff, ctcf_out, xlo, xhi, ylo, yhi, zlo, zhi):

    random.seed(seed)

    bead_stats, unique_ids = sequence_info(sequence, ctcf)
    n_molecules = sum(bead_stats.values())
    type_list = list(bead_stats.keys())
    n_types = len(type_list)
    write_header(n_molecules, n_types, xlo, xhi, ylo, yhi, zlo, zhi)
    write_masses(type_list)
    write_atoms(sequence, unique_ids, type_list)
    write_bonds(n_molecules)
    write_angles(n_molecules)
    write_ctcf_interactions(type_list, ctcf_coeff, ctcf_out)


def write_ctcf_interactions(type_list, ctcf_coeff, ctcf_out):
    ctcf_coeff = ctcf_coeff.replace(',',' ')
    with pct.open(ctcf_out, mode='w', stderr=False) as f:
        for idx1, id1 in enumerate(type_list):
            if not id1.startswith('F'):
                continue
            for id2 in type_list[idx1 + 1:]:
                if not id2.startswith('R'):
                    continue
                idx2 = type_list.index(id2)
                f.write(f'pair_coeff {idx1+1} {idx2+1} {ctcf_coeff}\n')


def sequence_info(path, ctcf):
    bead_stats = defaultdict(int)
    unique_ids = []
    with pct.open(path) as f:
        for line in f:
            bead = line.strip()
            if ctcf and bead in ['F', 'R']:
                id = len(unique_ids) + 1
                bead = f'{bead}-{id}'
                unique_ids.append(bead)
            bead_stats[bead] += 1
    return bead_stats, unique_ids


def write_header(n_molecules, n_types, xlo, xhi, ylo, yhi, zlo, zhi):
    sys.stdout.write(
        'LAMMPS data file from restart file: timestep = 0, procs = 1\n\n'
        f'{n_molecules} atoms\n'
        f'{n_molecules-1} bonds\n'
        f'{n_molecules-2} angles\n\n'
        f'{n_types} atom types\n'
        f'1 bond types\n'
        f'1 angle types\n\n'
        f'{xlo} {xhi} xlo xhi\n'
        f'{ylo} {yhi} ylo yhi\n'
        f'{zlo} {zhi} zlo zhi\n\n')


def write_masses(type_list):
    sys.stdout.write('Masses\n\n')
    for n, type in enumerate(type_list):
        sys.stdout.write(f'{n+1} 1 # {type}\n')
    sys.stdout.write('\n')


def write_atoms(sequence, unique_ids, type_list, r=1.1):
    sys.stdout.write('Atoms\n\n')
    with pct.open(sequence) as f:
        for i, line in enumerate(f):
            type = line.strip()
            if unique_ids and type in ['F', 'R']:
                type = unique_ids[0] # Get next first unique_id
                unique_ids = unique_ids[1:] # Remove from list
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
            atom_number = type_list.index(type) + 1
            sys.stdout.write(f'{i+1} 1 {atom_number} {x} {y} {z} 0 0 0\n')
            prev_x = x
            prev_y = y
            prev_z = z
    sys.stdout.write('\n')


def write_bonds(n_molecules):
    sys.stdout.write('Bonds\n\n')
    n_bonds = n_molecules - 1
    for b in range(1, n_bonds + 1):
        sys.stdout.write(f'{b} 1 {b} {(b%n_molecules)+1}\n')
    sys.stdout.write('\n')


def write_angles(n_molecules):
    sys.stdout.write('Angles\n\n')
    n_angles = n_molecules - 2
    for a in range(1, n_angles + 1):
        sys.stdout.write(f'{a} 1 {a} {a+1} {a+2}\n')
    sys.stdout.write('\n')


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
