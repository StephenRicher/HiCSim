#!/usr/bin/env python3

import sys
import json
import random
import logging
import argparse
import pandas as pd
from mpi4py import MPI
from lammps import lammps
from utilities import readJSON
from collections import namedtuple, defaultdict
from argUtils import setDefaults, createMainParent


__version__ = '1.0.0'


def runLammps(equil: str, atomGroups: str, simTime: int, TADStatus: str,
              sepThresh: float, nExtrudersPerBead: float, extrusionRate: float,
              updateInterval: float, onRate: float, offRate: float,
              writeInterval: float, extrusion: bool, harmonicCoeff: float,
              timestep: float, TFswap: float, radiusGyrationOut: str,
              pairCoeffs: str, simOut: str, beadTypes: str, nBasesPerBead: int,
              seed: int):

    random.seed(seed)

    offRate = computeOffRate(nBasesPerBead, extrusionRate, prob=0.5, distance=40_000)
    stepProb = updateInterval * extrusionRate # Advancement probability
    addProb  = updateInterval * onRate        # Attachment probability
    offProb  = updateInterval * offRate       # Dettachment probability

    # Convert time unit to timesteps
    writeIntervalSteps = int(writeInterval / timestep)
    TFswapSteps = int(TFswap / timestep)

    global lmp
    lmp = lammps()

    lmp.command(f'read_restart {equil}')
    lmp.command('reset_timestep 0')

    lmp.command('neigh_modify every 1 delay 1 check yes')
    lmp.command('fix 1 all nve')
    lmp.command(f'fix 2 all langevin 1.0 1.0 1.0 {seed}')

    atomGroups = readJSON(atomGroups)

    # Assign all groups
    writeGroups(atomGroups)
    # Assign atom types
    beadTypes = readJSON(beadTypes)
    for type, typeID in beadTypes.items():
        if (type in ['TFa', 'TFi', 'N']) or (type not in atomGroups):
            continue
        lmp.command(f'set group {type} type {typeID}')

    # Reduce timestep and run short equilibration before setting pair coeffs
    lmp.command(f'timestep {timestep / 10}')

    lmp.command('run 100000') # step 1

    lmp.command('bond_style hybrid fene harmonic')
    lmp.command('special_bonds fene')
    lmp.command('bond_coeff 1 fene 30 1.5 1.0 1.0')
    lmp.command('bond_coeff 2 fene 30 7.0 1.0 1.0')
    #lmp.command(f'bond_coeff 2 harmonic {harmonicCoeff} 1.5')

    lmp.command('pair_style lj/cut 1.12246')
    lmp.command('pair_modify shift yes')
    lmp.command('pair_coeff * *  1.0 1.0 1.12246') # Repulsive - everything with everything
    
    lmp.command('run 100000') # step 2

    # Increase timestep back
    lmp.command(f'timestep {timestep}')

    if pairCoeffs:
        writePairCoeffs(pairCoeffs, beadTypes)

    lmp.command('compute RG DNA gyration')
    lmp.command('variable RG equal c_RG')
    lmp.command(f'fix RG DNA print {writeIntervalSteps} "$(step) ${{RG}}" file {radiusGyrationOut}')

    lmp.command(f'dump 2 all custom {writeIntervalSteps} {simOut} id type x y z ix iy iz')
    lmp.command("dump_modify 2 format line '%d %d %.5f %.5f %.5f %d %d %d' sort 1")

    lmp.command('reset_timestep 0')
    lmp.command('neighbor 5 bin')

    lmp.command(f'fix swap TU atom/swap {TFswapSteps} 5 {seed} 10 ke yes types 1 2')

    allStatus = []
    allExtruders = Extruders(
        nExtrudersPerBead, atomGroups, offProb, stepProb, addProb, sepThresh)

    nIntervals = int(simTime / updateInterval)
    updateIntervalSteps = int(updateInterval / timestep)
    for step in range(nIntervals):
        status = allExtruders.writeTADs()
        status['timestep'] = step * updateIntervalSteps
        allStatus.append(status)
        if step > 5:
            allExtruders.updateExtrusion()
        lmp.command(f'run {updateIntervalSteps}')
    lmp.command('undump 2')
    lmp.command('unfix RG')
    # Update TAD status for last timestep
    status = allExtruders.writeTADs()
    status['timestep'] = (step + 1) * updateIntervalSteps
    allStatus.append(status)
    allStatus = pd.concat(allStatus)
    allStatus['time'] = allStatus['timestep'] * timestep
    allStatus.to_csv(TADStatus, header=True, index=False)


def computeOffRate(nBasesPerBead, extrusionRate, prob=0.5, distance=40_000):
    """ Compute per time offProb required such that extruded
        distance is reached when the cumulative offProb reaches prob """
    # Number of extruded beads to reach extrusion distance
    extrusionBeads = distance / nBasesPerBead
    # Number steps required to reach that distance, given stepProb
    totalSteps = extrusionBeads / extrusionRate
    # Compute offProb such that there is a 50% offRate after that step number
    offProb = 1 - (prob ** (1 / totalSteps))
    return offProb


def writePairCoeffs(pairCoeffs, beadType):
    """ Define typeID pairCoeffs from beadType mappings"""
    with open(pairCoeffs) as fh:
        for line in fh:
            line = line.strip().split()
            coeff = ' '.join(line[2:])
            type1 = line[0]
            type2 = line[1]
            try:
                type1ID = '*' if type1 == '*' else beadType[type1]
                type2ID = '*' if type2 == '*' else beadType[type2]
            except KeyError:
                logging.error(
                    f'Atleast one of {type1} or {type2} does not exist. Skipping.')
                continue
            if (type1ID != '*') and (type2ID != '*'):
                # Ensure the smaller typeID is written first
                if type1ID > type2ID:
                    temp = type1ID
                    type1ID = type2ID
                    type2ID = temp
            lmp.command(f'pair_coeff {type1ID} {type2ID} {coeff}\n')


def writeGroups(atomGroups):
    """ Read lammps input and define bead ID which are polymer """
    for group, beadIDs in atomGroups.items():
        if (group == 'DNA') or (group == 'TF'):
            continue # Already defined in equilibration
        ids = " ".join([str(i) for i in beadIDs])
        lmp.command(f'group {group} id {ids}')


class Extruder():
    def __init__(self):
        self.isBound = False
        self.left = None
        self.right = None

    @property
    def left(self):
        return self._left

    @left.setter
    def left(self, left):
        self._left = left

    @property
    def right(self):
        return self._right

    @right.setter
    def right(self, right):
        self._right = right

    def loopIDs(self):
        return range(self.left, self.right + 1)

    def nextPosition(self, direction, N=1):
        assert direction in ['left', 'right']
        if direction == 'left':
            return self.currentPosition(direction) - N
        else:
            return self.currentPosition(direction) + N


    def currentPosition(self, direction):
        assert direction in ['left', 'right']
        if not self.isBound:
            return None
        elif direction == 'left':
            return self.left
        else:
            return self.right


    def attach(self, left, right, sepThresh):
        if self.isBound:
            logging.warning(f'Extruding from {self.left}:{self.right} to {left}:{right}')
            # Delete old bonds before creating new.
            self.detach(message=False)
        else:
            logging.warning(f'Attaching to {left}:{right}')
        self.isBound = True
        self.left = left
        self.right = right
        # Create harmonic between new bead pairs
        name = f'{self.left}-{self.right}'
        lmp.command(f'group {name} id {self.left} {self.right}')
        lmp.command(f'create_bonds many {name} {name} 2 0.0 {sepThresh}')


    def detach(self, message=True):
        name = f'{self.left}-{self.right}'
        lmp.command(f'delete_bonds {name} bond 2 remove special')
        lmp.command(f'group {name} delete')
        if message:
            logging.warning(f'Detaching {self.left}:{self.right}')
        self.isBound = False
        self.left = None
        self.right = None


    def occupied(self):
        return [self.left, self.right]


class Extruders():
    def __init__(self, nExtrudersPerBead, atomGroups, removeProb, stepProb, addProb, sepThresh):
        self.beadIDs = defaultdict(list, atomGroups)
        self.nExtruders = int(len(self.beadIDs['DNA']) * nExtrudersPerBead)
        self.extruders = [Extruder() for i in range(self.nExtruders)]
        self.removeProb = removeProb
        self.stepProb = stepProb
        self.addProb = addProb
        self.sepThresh = sepThresh
        self.coords = None
        self.box = self.processBox(lmp.extract_box())
        self.TADstatus = None
        # Variables to track whether the previous status was updated
        self.removed = False
        self.updated = False
        self.updated = False


    def processBox(self, lmpBox):
        """ Convert output of extract_box() to box dimensions """
        box = namedtuple('box', 'xSize ySize zSize')
        xSize = abs(lmpBox[0][0] - lmpBox[1][0])
        ySize = abs(lmpBox[0][1] - lmpBox[1][1])
        zSize = abs(lmpBox[0][2] - lmpBox[1][2])
        return box(xSize, ySize, zSize)


    def isForwardCTCF(self, beadID):
        """ Retrun True if beadID is a 'F' or 'B' CTCF site """
        return (beadID in self.beadIDs['B']) or (beadID in self.beadIDs['F'])


    def isReverseCTCF(self, beadID):
        """ Retrun True if beadID is a 'R' or 'B' CTCF site """
        return (beadID in self.beadIDs['B']) or (beadID in self.beadIDs['R'])


    def isCTCF(self, beadID):
        """ Return True if beadID is CTCF """
        return self.isForwardCTCF(beadID) or self.isReverseCTCF(beadID)


    def farBeadID(self, direction):
        assert direction in ['left', 'right']
        if direction == 'left':
            return min(self.beadIDs['DNA'])
        else:
            return max(self.beadIDs['DNA'])


    def outOfRange(self, beadID):
        """ Return True if beadID is outside of DNA beadID range """
        return (beadID < self.farBeadID('left')) or (beadID > self.farBeadID('right'))


    def inLoop(self, pos):
        """ Check if beadID is within an extrusion loop """
        nLoops = 0
        for extruder in self.boundExtruders():
            if pos in extruder.loopIDs():
                nLoops += 1
        return nLoops


    def isOccupied(self, beadID):
        """ Check if beadID corresponds to an occupied bead """
        for extruder in self.boundExtruders():
            if beadID in extruder.occupied():
                return True
        return False


    def isValid(self, beadID):
        """ Return True if beadID is a valid cohesin binding position """
        if (self.isOccupied(beadID)
                or self.isCTCF(beadID)
                or self.outOfRange(beadID)):
            return False
        else:
            return True


    def boundExtruders(self):
        allBound = []
        for extruder in self.extruders:
            if extruder.isBound:
                allBound.append(extruder)
        return allBound


    def unboundExtruders(self):
        allUnbound = []
        for extruder in self.extruders:
            if not extruder.isBound:
                allUnbound.append(extruder)
        return allUnbound


    def removeExtruders(self):
        rc = False
        for extruder in self.boundExtruders():
            if random.random() > self.removeProb:
                continue
            extruder.detach()
            rc = True
        return rc


    def canExtrude(self, extruder, direction):
        """ Check if extruder can validly extrude in direction """
        assert direction in ['left', 'right']
        # Retrieve next position1 in relevant direction
        nextPos = extruder.nextPosition(direction, 1)
        if self.isOccupied(nextPos):
            return False
        currentPos = extruder.currentPosition(direction)
        if direction == 'left' and self.isForwardCTCF(currentPos):
            return False
        if direction == 'right' and self.isReverseCTCF(currentPos):
            return False
        if random.random() > self.stepProb:
            return False
        return True


    def updateExtruders(self):
        rc = False
        for extruder in self.boundExtruders():
            next = {}
            updated = False
            for direction in ['left', 'right']:
                if self.canExtrude(extruder, direction):
                    next[direction] = extruder.nextPosition(direction)
                    updated = True
                else:
                    next[direction] = extruder.currentPosition(direction)
            # Skip if no extrusion
            if not updated:
                continue
            # Detach extruder if the next position slides off edge of polymer
            if self.outOfRange(next['left']) or self.outOfRange(next['right']):
                extruder.detach()
                continue
            if self.beadSeperation(next['left'], next['right']) > self.sepThresh:
                continue
            extruder.attach(next["left"], next["right"], self.sepThresh)
            rc = True
        return rc


    def attachExtruders(self):
        rc = False
        for extruder in self.unboundExtruders():
            if random.random() > self.addProb:
                continue
            startPos = random.choice(self.beadIDs['DNA'])
            endPos = startPos + 2
            # Check position before and position after strong bond
            for pos in range(startPos, endPos + 1):
                if not self.isValid(pos):
                    break
            else:
                if self.beadSeperation(startPos, endPos) > self.sepThresh:
                    continue
                extruder.attach(startPos, endPos, self.sepThresh)
                rc = True
        return rc


    def updateExtrusion(self):
        """ Perform extrusion and store whether changes were made """
        self.removed = self.removeExtruders()
        self.updated = self.updateExtruders()
        self.attached = self.attachExtruders()


    def anyUpdated(self):
        """ Return True if any extruders were removed/updated or attached """
        return self.removed or self.updated or self.attached


    def updateCoordinates(self):
        """ Store all Atom coordinates keyed by atom ID """
        pos = namedtuple('pos', 'x y z')
        allPos = {}
        coords = iter(lmp.gather_atoms("x", 1, 3))
        for ID, (x, y, z) in enumerate(zip(coords, coords, coords), 1):
            allPos[ID] = pos(x, y, z)
        self.coords = allPos


    def beadSeperation(self, beadID1, beadID2):
        self.updateCoordinates()
        beadID1pos = self.coords[beadID1]
        beadID2pos = self.coords[beadID2]
        total = 0
        for i, (a, b) in enumerate(zip(beadID1pos, beadID2pos)):
            delta = abs(b - a)
            if delta > self.box[i] - delta:
                delta = self.box[i] - delta
            total += delta ** 2
        return total ** 0.5


    def writeTADs(self):
        """ Write per-bead TAD status and section identifier
            Section ID is unique per timepoint and groups beads
            within a pair of extruder boundaries """
        if (self.TADstatus is None) or (self.anyUpdated):
            allStatus = {'ID': [], 'TADstatus': [], 'TADgroup': []}
            TADgroup = 0 # Track number of boundaries crossed
            for bead in self.beadIDs['DNA']:
                if self.isOccupied(bead):
                    status = -1
                    TADgroup += 1
                else:
                    # Number of extrusion loops the bead overlaps
                    # Can be > 1 if an extruder binds to an extruded loop
                    status = self.inLoop(bead)
                allStatus['ID'].append(bead)
                allStatus['TADstatus'].append(status)
                allStatus['TADgroup'].append(TADgroup)
            self.TADstatus = allStatus
        return pd.DataFrame(self.TADstatus)


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=runLammps)
    parser.add_argument('equil', help='Lammps equilibration restart file.')
    parser.add_argument('atomGroups', help='Atom groups JSON file from polymer')
    parser.add_argument(
        '--simTime', type=int, default=100000,
        help='Total simulation time (default: %(default)s)')
    parser.add_argument(
        '--sepThresh', type=float, default=6.0,
        help='Max distance for bond creation (default: %(default)s)')
    parser.add_argument(
        '--nExtrudersPerBead', type=float, default=(30.0 / 5000.0),
        help='Number of extruders per polymer bead (default: %(default)s)')
    parser.add_argument(
        '--extrusionRate', type=float, default=(1 / 500),
        help='Time per extrusion, higher interval = slower extrusion '
             '(default: %(default)s)')
    parser.add_argument(
        '--updateInterval', type=float, default=10.0,
        help='Interval between extrusion status update (default: %(default)s)')
    parser.add_argument(
        '--onRate', type=float, default=(1.0 / 1000.0),
        help='Rate at which extruders attach in inverse timesteps '
             '(default: %(default)s)')
    parser.add_argument(
        '--offRate', type=float, default=(1.0 / 40_000.0),
        help='Rate at which extruders detach in inverse timesteps '
             '(default: %(default)s)')
    parser.add_argument(
        '--TADStatus', default=sys.stdout,
        help='Path to write TAD status of each bead per '
             'update-time (default: stdout)')
    parser.add_argument(
        '--writeInterval', type=float, default=10,
        help='Simulation status output interval in time units '
             '(default: %(default)s)')
    parser.add_argument(
        '--extrusion', default=False, action='store_true',
        help='Set to prepare script for loop extrusion (default: %(default)s)')
    parser.add_argument(
        '--harmonicCoeff', type=float, default=40.0,
        help='Harmonic bond strength for extrusion factors '
             '(default: %(default)s)')
    parser.add_argument(
        '--TFswap', type=float, default=100,
        help='Swap frequency of TF active/inactivation in time units '
             '(default: %(default)s)')
    parser.add_argument(
        '--seed', type=int, default=random.randint(0, 10e9),
        help='Seed for simulation (default: %(default)s)')
    parser.add_argument(
        '--pairCoeffs', help='Bead pair coefficient definitions.')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '--radiusGyrationOut', required=True,
        help='Output file for radius of gyration data.')
    requiredNamed.add_argument(
        '--simOut', required=True, help='Output file for simulation.')
    requiredNamed.add_argument(
        '--timestep', type=float, required=True,
        help='Size of timestep in time units from equil.')
    requiredNamed.add_argument(
        '--beadTypes', required=True,
        help='Mapping of bead type to type ID from equil.')
    requiredNamed.add_argument(
        '--nBasesPerBead', type=int, required=True,
        help='Number of bases that represent a single bead.')
    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
