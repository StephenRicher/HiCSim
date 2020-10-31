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


class Sequence:

    def __init__(self, sequence, basesPerBead, atac=None,
            name='DNA', ctcf=False, beadID=1):
        self._beadID = beadID
        self.basesPerBead = basesPerBead
        # Read ATAC modifiers for each Bead
        if atac is not None:
            self.atac = readJSON(atac)
        else:
            self.atac = defaultdict(int)
        # Dictionary to store bead type : modifier
        self.modifiers = {}
        # Should beads set as 'F' or 'R' be processed as CTCF?
        self._ctcf = ctcf
        self.ctcfs = {}
        self.name = name
        self.load(sequence)

    def load(self, sequence):
        self.sequence = {}
        self.types = []
        with open(sequence) as f:
            for i, line in enumerate(f):
                type = line.strip()
                if self._ctcf and type in ['F', 'R', 'B']:
                    type = f'{type}-{self._beadID}'
                    # Associate unique CTCF type with bead ID
                    self.ctcfs[type] = self._beadID
                if type == 'N':
                    type = f'{type}-{self._beadID}'
                    self.modifiers[type] = float(self.atac[str(i)])
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

    def __init__(self, nMonomers=0, polymerSeed=random.randint(1, 1e100),
            monomerSeed=random.randint(1, 1e100), randomWalk=False, TU_types=[]):
        self.sequences = []
        self.TU_types = TU_types # Beads to be considered transcriptional units
        self._beadID = 1
        # typeID 1 reserved for active TFs
        # typeID 2 for inactive TFs
        self._typeID = 3
        self._bondID = 1
        self.typeIDs = {}
        self.randomWalk = randomWalk
        self._monomerSeed = monomerSeed
        self._polymerSeed = polymerSeed
        self.monomers = self.loadMonomer(nMonomers, prop=0.5)


    def loadSequence(self, sequence_file, atac, basesPerBead, ctcf):
        sequence = Sequence(sequence_file, basesPerBead=basesPerBead, atac=atac,
            ctcf=ctcf, beadID=self._beadID)
        self.sequences.append(sequence)
        self._beadID += sequence.nBeads
        for type in sequence.types:
            if type in ['TFa', 'TFi']:
                sys.exit('Beads "TFa" and "TFi" reserved for monomer beads.')
            elif type == 'None':
                sys.exit('Cannot set Bead type as "None".')
            self.addType(type)


    def loadMonomer(self, nMonomers, prop=0.5):
        # Set monomer number
        self.nMonomers = nMonomers
        # Prop is proportion of active monomers
        # List of active monomers ('TFa') and inactive monomers ('TFi')
        monomers = ( int(nMonomers * prop)       * ['TFa']
                   + int(nMonomers * (1 - prop)) * ['TFi'])
        monomerIDs = {}
        for monomer in monomers:
            self.addType(monomer)
            monomerIDs[self._beadID] = monomer
            self._beadID += 1
        return monomerIDs


    def loadBox(self, xlo, xhi, ylo, yhi, zlo, zhi):
        Dim = namedtuple('Dim', 'lo hi')
        Box = namedtuple('Box', ['x', 'y', 'z'])
        self.box = Box(Dim(xlo, xhi), Dim(ylo, yhi), Dim(zlo, zhi))


    def addType(self, type):
        # Only reserve types 1 and 2 if monomers are in simulation ...
        # Bead 'TFa' is reserved for active monomrs - type 1
        if self.nMonomers and type == 'TFa':
            self.typeIDs[type] = 1
        # Bead 'TFi' is reserved for inactive monomrs - type 2
        elif self.nMonomers and type == 'TFi':
            self.typeIDs[type] = 2
        elif type not in self.typeIDs:
            self.typeIDs[type] = self._typeID
            self._typeID += 1

    @property
    def nBeads(self):
        """ Return total number of beads in simulated. """
        n = sum([sequence.nBeads for sequence in self.sequences])
        n += self.nMonomers
        return n


    @property
    def nAngles(self):
        return sum([sequence.nAngles for sequence in self.sequences])

    @property
    def nBonds(self):
        bonds = 0
        for sequence in self.sequences:
            bonds += sequence.nBonds + len(self.detectConvertCTCF(sequence))
        return bonds

    @property
    def nBondTypes(self):
        for sequence in self.sequences:
            if len(self.detectConvertCTCF(sequence)) > 0:
                n = 2
                break
        else:
            n = 1
        return n

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
            f'{self.nBondTypes} bond types\n'
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
        for sequence in self.sequences:
            if self.randomWalk:
                x, y, z = self.makeRandomWalk(sequence.nBeads)
            else:
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
                TU = 'TU' if type in self.TU_types else ''
                sys.stdout.write(f'{beadID} 1 {typeID} {x[i]} {y[i]} {z[i]} 0 0 0 # DNA {type} {TU}\n')


    def writeMonomers(self):
        random.seed(self._monomerSeed)
        for beadID, type in self.monomers.items():
            x = random.uniform(self.box.x.lo, self.box.x.hi)
            y = random.uniform(self.box.y.lo, self.box.y.hi)
            z = random.uniform(self.box.z.lo, self.box.z.hi)
            typeID = self.typeIDs[type]
            sys.stdout.write(f'{beadID} 1 {typeID} {x} {y} {z} 0 0 0 # TF\n')


    def writeBonds(self):
        sys.stdout.write('Bonds\n\n')
        for sequence in self.sequences:
            nBonds = sequence.nBonds
            for i, beadID in enumerate(sequence.sequence):
                if i == sequence.nBonds:
                    break
                sys.stdout.write(f'{self._bondID} 1 {beadID} {beadID+1}\n')
                self._bondID += 1
            for forward, reverse in self.detectConvertCTCF(sequence):
                sys.stdout.write(f'{self._bondID} 2 {forward} {reverse}\n')
                self._bondID += 1
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
        typeIdxs = list(range(len(typeList) - 1))
        # Ensure shuffle is same for each call within object
        random.seed(self._polymerSeed)
        # Shuffle to initiate loop extrusion differently per replicate
        random.shuffle(typeIdxs)
        usedSites = []
        for typeIdx in typeIdxs:
            # Set to True to prevent loops forming passed already formed loops
            alreadyUsed = False
            # Reverse down the list until reach forward CTCF
            for forwardBead in reversed(typeList[:typeIdx + 1]):
                if forwardBead in usedSites:
                    alreadyUsed = True
                if forwardBead.startswith(('F', 'B')) and not alreadyUsed:
                    break
            else:
                continue
            # Move up the list until reach reverse CTCF
            for reverseBead in typeList[typeIdx + 1:]:
                if reverseBead in usedSites:
                    alreadyUsed = True
                if reverseBead.startswith(('R', 'B')) and not alreadyUsed:
                    break
            else:
                continue
            usedSites.extend([forwardBead, reverseBead])
            pairs.add((sequence.ctcfs[forwardBead], sequence.ctcfs[reverseBead]))
        return pairs


    def writeCoeffs(self, type1, type2, coeff, fh=sys.stderr):
        try:
            type1ID = '*' if type1 == '*' else self.typeIDs[type1]
            type2ID = '*' if type2 == '*' else self.typeIDs[type2]
        except KeyError:
            logging.error(
                f'Atleast one of {type1} or {type2} does not exist. Skipping.')
            return 1
        if (type1ID != "*") and (type2ID != "*"):
            # Ensure the smaller typeID is written first
            if type1ID > type2ID:
                temp = type1ID
                type1ID = type2ID
                type2ID = temp
        fh.write(f'pair_coeff {type1ID} {type2ID} {coeff}\n')


    def writeATACCoeffs(self, TFcoeffMin, TFcoeffMax, fh=sys.stderr):
        """ Write ATAC modified pair coefficients for N beads and monomers """

        # Range between min and max coefficients
        coeffRange = TFcoeffMax - TFcoeffMin
        TFa_typeID = 1
        for sequence in self.sequences:
            for type, modifier in sequence.modifiers.items():
                N_typeID = self.typeIDs[type]
                TFcoeffMod = (modifier * coeffRange) + TFcoeffMin
                fh.write(
                    f'pair_coeff {TFa_typeID} {N_typeID} {TFcoeffMod} 1.0 1.8\n')



    def writeGroup(self, fh=sys.stderr):
        for sequence in self.sequences:
            firstBeadID = min(sequence.sequence)
            lastBeadID = max(sequence.sequence)
            name = sequence.name
            fh.write(f'group {name} id {firstBeadID}:{lastBeadID}\n')
        fh.write(f'group monomer id {min(self.monomers)}:{max(self.monomers)}\n')
        #fh.write(f'group monomer id {lastBeadID + 1}:{lastBeadID + self.nMonomers}\n')



def main(file, atac, monomerSeed, polymerSeed, nMonomers, ctcf, basesPerBead, randomWalk,
        pairCoeffs, coeffOut, groupOut, xlo, xhi, ylo, yhi, zlo, zhi, **kwargs):

    # TO DO - the tye list needs to be input by user
    TU_types=['1', '7', '3']
    dat = lammps(nMonomers, monomerSeed=monomerSeed, polymerSeed=polymerSeed,
        randomWalk=randomWalk, TU_types=TU_types)
    dat.loadBox(xlo, xhi, ylo, yhi, zlo, zhi)
    dat.loadSequence(file, atac, basesPerBead, ctcf)
    dat.writeLammps()
    with open(pairCoeffs) as fh:
        with open(coeffOut, 'w') as out:
            for line in fh:
                line = line.strip().split()
                coeff = ' '.join(line[2:])
                atom1 = line[0]
                atom2 = line[1]
                # IF F-R or R-C (CTCF interaction)
                dat.writeCoeffs(atom1, atom2, coeff, fh=out)
            dat.writeATACCoeffs(TFcoeffMin=3,  TFcoeffMax=8, fh=out)
    with open(groupOut, 'w') as out:
        dat.writeGroup(fh=out)


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
        'file',metavar='SEQUENCE', help='Input bead sequence file.')
    custom.add_argument(
        '--atac', help='ATAC bead binding modifier')
    custom.add_argument(
        '--nMonomers', default=0, type=int,
        help='Set number of monomers for simulation.')
    custom.add_argument(
        '--ctcf', default=False, action='store_true',
        help='Process beads labelled F and R as CTCF sites. Each '
             'bead is given a unique ID and convergent orientation '
             'pair coeffs are output (default: %(default)s)')
    custom.add_argument(
        '--pairCoeffs',
        help='Input of bead name pair coefficients.')
    custom.add_argument(
        '--coeffOut', default=None,
        help='Outfile to write atom ID pair coefficients (default: stderr)')
    custom.add_argument(
        '--groupOut', default=None,
        help='Outfile to write atom group ID assignments (default: stderr)')
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
