#!/usr/bin/env python3

import sys
import json
import random
import logging
import argparse
import pandas as pd
from lammps import lammps
from collections import namedtuple, defaultdict
from argUtils import setDefaults, createMainParent


__version__ = '1.0.0'


def runLammps(initLam: str, atomGroups: str, simTime: int, TADStatus: str,
              sepThresh: float, nExtrudersPerBead: float, extrusionRate: float,
              updateInterval: float, onRate: float, offRate: float):

    global lmp
    lmp = lammps()
    lmp.file(initLam)

    dt = getTimestep(initLam)
    nIntervals = int(simTime / updateInterval)
    stepProb = updateInterval * extrusionRate # Advancement probability
    addProb  = updateInterval * onRate        # Attachment probability
    offProb  = updateInterval * offRate       # Dettachment probability

    allStatus = {}
    allExtruders = Extruders(
        nExtrudersPerBead, atomGroups, offProb, stepProb, addProb, sepThresh)
    for step in range(nIntervals):
        if step > 0:
            allExtruders.updateExtrusion()
        time = step * updateInterval
        allStatus[time] = allExtruders.writeTADs()
        lmp.command(f'run {int(updateInterval / dt)}')

    allStatus = (pd.DataFrame(
        allStatus, index=allExtruders.beadIDs['DNA'])
        .reset_index().melt(id_vars='index'))
    allStatus.columns = ['beadID', 'time', 'tadStatus']
    allStatus.to_csv(TADStatus, header=True, index=False)


def getTimestep(file):
    """ Retrieve timestep from lammps init file """
    with open(file) as fh:
        for line in fh:
            line = line.strip().split()
            if (len(line) == 2) and (line[0] == 'timestep'):
                return float(line[1])


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

    def nextPosition(self, direction):
        assert direction in ['left', 'right']
        move = {'left' : -1, 'right' : 1}
        return self.currentPosition(direction) + move[direction]

    def currentPosition(self, direction):
        assert direction in ['left', 'right']
        if not self.isBound:
            return None
        elif direction == 'left':
            return self.left
        else:
            return self.right

    def attach(self, left, right):
        self.isBound = True
        self.left = left
        self.right = right

    def detach(self):
        self.isBound = False
        self.left = None
        self.right = None


class Extruders():
    def __init__(self, nExtrudersPerBead, beadIDs, removeProb, stepProb, addProb, sepThresh):
        self.beadIDs = self.processAtomGroup(beadIDs)
        self.nExtruders = int(len(self.beadIDs['DNA']) * nExtrudersPerBead)
        self.extruders = [Extruder() for i in range(self.nExtruders)]
        self.removeProb = removeProb
        self.stepProb = stepProb
        self.addProb = addProb
        self.sepThresh = sepThresh
        self.coords = None
        self.box = self.processBox(lmp.extract_box())


    def processBox(self, lmpBox):
        """ Convert output of extract_box() to box dimensions """
        box = namedtuple('box', 'xSize ySize zSize')
        xSize = abs(lmpBox[0][0] - lmpBox[1][0])
        ySize = abs(lmpBox[0][1] - lmpBox[1][1])
        zSize = abs(lmpBox[0][2] - lmpBox[1][2])
        return box(xSize, ySize, zSize)


    def processAtomGroup(self, beadIDs):
        """ Process JSON atom groups to store CTCF and DNA bead IDs """
        with open(beadIDs) as fh:
            return defaultdict(list, json.load(fh))


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
            if beadID in [extruder.left, extruder.right]:
                return True
        else:
            return False


    def isValid(self, beadID):
        """ Return True if beadID is a valid cohesin binding position """
        if self.isOccupied(beadID) or self.isCTCF(beadID) or self.outOfRange(beadID):
            return False
        else:
            return True


    def canExtrude(self, extruder, direction):
        """ Check if extruder can validly extrude in direction """
        assert direction in ['left', 'right']
        # Check if next position is already occupied
        if self.isOccupied(extruder.nextPosition(direction)):
            return False
        currentPos = extruder.currentPosition(direction)
        if direction == 'left' and self.isForwardCTCF(currentPos):
            return False
        if direction == 'right' and self.isReverseCTCF(currentPos):
            return False
        if random.random() > self.stepProb:
            return False
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
        for extruder in self.boundExtruders():
            if random.random() > self.removeProb:
                continue
            lmp.command(f'group newbond id {extruder.left} {extruder.right}')
            lmp.command('delete_bonds newbond bond 2 remove special')
            lmp.command('group newbond delete')
            extruder.detach()


    def updateExtruders(self):
        for extruder in self.boundExtruders():
            next = {}
            for direction in ['left', 'right']:
                updated = False
                if self.canExtrude(extruder, direction):
                    next[direction] = extruder.nextPosition(direction)
                    updated = True
                else:
                    next[direction] = extruder.currentPosition(direction)
            # Skip if no extrusion
            if not updated:
                continue
            if self.beadSeperation(next['left'], next['right']) > self.sepThresh:
                continue
            lmp.command(f'group tmpUpdate id {extruder.left} {extruder.right}')
            lmp.command('delete_bonds tmpUpdate bond 2 remove special')
            lmp.command('group tmpUpdate delete')
            # Detach extruder if the next position slides of edge of polymer
            if self.outOfRange(next['left']) or self.outOfRange(next['right']):
                extruder.detach()
                continue
            lmp.command(f'group tmpUpdate id {next["left"]} {next["right"]}')
            lmp.command(f'create_bonds many tmpUpdate tmpUpdate 2 0.0 {self.sepThresh}')
            lmp.command('group tmpUpdate delete')
            extruder.left = next['left']
            extruder.right = next['right']


    def attachExtruders(self):
        for extruder in self.unboundExtruders():
            if random.random() > self.addProb:
                continue
            startPos = random.choice(self.beadIDs['DNA'])
            for pos in range(startPos, startPos + 2 + 1):
                if not self.isValid(pos):
                    break
            else:
                if self.beadSeperation(startPos, startPos + 2) < self.sepThresh:
                    extruder.attach(startPos, startPos + 2)
                    lmp.command(f'group newbond id {extruder.left} {extruder.right}')
                    lmp.command(f'create_bonds many newbond newbond 2 0.0 {self.sepThresh}')
                    lmp.command('group newbond delete')


    def updateExtrusion(self):
        self.removeExtruders()
        self.updateExtruders()
        self.attachExtruders()


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
        """ Write per-bead TAD status """
        allStatus = []
        for bead in self.beadIDs['DNA']:
            if self.isOccupied(bead):
                status = -1
            else:
                status = self.inLoop(bead)
            allStatus.append(status)
        return allStatus


def parseArgs():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    mainParent = createMainParent(verbose=False, version=__version__)
    parser = argparse.ArgumentParser(
        epilog=epilog, description=__doc__, parents=[mainParent])
    parser.set_defaults(function=runLammps)
    parser.add_argument('initLam', help='Lammps initiation script.')
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
        '--offRate', type=float, default=(1.0 / 40000.0),
        help='Rate at which extruders detach in inverse timesteps '
             '(default: %(default)s)')
    parser.add_argument(
        '--TADStatus', default=sys.stdout,
        help='Path to write TAD status of each bead per '
             'update-time (default: stdout)')

    return setDefaults(parser)


if __name__ == '__main__':
    args, function = parseArgs()
    sys.exit(function(**vars(args)))
