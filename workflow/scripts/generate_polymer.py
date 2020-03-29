#!/usr/bin/env python3

""" Script to generate a random linear polymer as an input file for LAMMPS """

import sys
import math
import random
import argparse
import pyCommonTools as pct

def main():

    __version__ = '1.0.0'

    parser = pct.make_parser(
        verbose=True, version=__version__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--n_molecules', default=1000, type=int,
        help='Number of molecules in polymer.')
    parser.add_argument(
        '--polymer_type', default='homopolymer',
        choices=['homopolymer', 'random_copolymer', 'clustered_copolymer'],
        help='Type of copolymer to generate.')
    parser.add_argument(
        '--xlo', default=-50, type=float,
        help='Lower x-axis of simuluation box.')
    parser.add_argument(
        '--xhi', default=50, type=float,
        help='Upper x-axis of simuluation box.')
    parser.add_argument(
        '--ylo', default=-50, type=float,
        help='Lower y-axis of simuluation box.')
    parser.add_argument(
        '--yhi', default=50, type=float,
        help='Upper y-axis of simuluation box.')
    parser.add_argument(
        '--zlo', default=-50, type=float,
        help='Lower z-axis of simuluation box.')
    parser.add_argument(
        '--zhi', default=50, type=float,
        help='Upper z-axis of simuluation box.')
    parser.add_argument(
        '--atoms', metavar="ATOM=MASS", nargs="+", action=ParseDict,
        default={'1':'1', '2':'1', '3':'1', '4':'1'},
        help='Provide atom type and mass as key=value pairs.')
    parser.set_defaults(function=generate_linear_polymer)

    return (pct.execute(parser))

def generate_linear_polymer(
        n_molecules, polymer_type, xlo, xhi, ylo, yhi, zlo, zhi, atoms):

    sys.stdout.write(
        'LAMMPS data file from restart file: timestep = 0, procs = 1\n\n'
        f'{n_molecules} atoms\n'
        f'{n_molecules-1} bonds\n'
        f'{n_molecules-2} angles\n\n'
        f'{len(atoms)} atom types\n'
        f'1 bond types\n'
        f'1 angle types\n\n'
        f'{xlo} {xhi} xlo xhi\n'
        f'{ylo} {yhi} ylo yhi\n'
        f'{zlo} {zhi} zlo zhi\n\n'
        f'Masses\n\n')

    for atom in atoms:
        sys.stdout.write(f'{atom} {atoms[atom]}\n')

    sys.stdout.write('\nAtoms\n\n')
    r = 1.1
    cluster_size = n_molecules / len(atoms)
    atom_type_idx = 0 # Set initial index for atom type dictionary
    type_boundary_idx = 0 # Set initial index of type boundary change
    for n in range(1, n_molecules+1):

        if polymer_type == 'homopolymer':
            type = list(atoms)[0]
        elif polymer_type == 'random_copolymer':
            type = random.choice(list(atoms))
        else:
            if n > type_boundary_idx + cluster_size:
                type_boundary_idx = n
                atom_type_idx += 1
            type = list(atoms)[atom_type_idx]

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

    sys.stdout.write('\nBonds\n\n')
    n_bonds = n_molecules - 1
    for b in range(1, n_bonds + 1):
        sys.stdout.write(f'{b} 1 {b} {(b%n_molecules)+1}\n')

    sys.stdout.write('\nAngles\n\n')
    n_angles = n_molecules - 2
    for a in range(1, n_angles + 1):
        sys.stdout.write(f'{a} 1 {a} {a+1} {a+2}\n')

class ParseDict(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        d = {}

        if values:
            for item in values:
                split_items = item.split("=", 1)
                key = split_items[0].strip()
                value = split_items[1]
                d[key] = value

        setattr(namespace, self.dest, d)

if __name__ == '__main__':
    sys.exit(main())
