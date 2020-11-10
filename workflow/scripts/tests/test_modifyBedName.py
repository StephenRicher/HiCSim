#!/usr/bin/env python

import os
import sys
import pytest
from testUtils import datadir, validateOutput
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from modifyBedName import *


modifyName_params = ([('test1.in', 'test1.out')])
@pytest.mark.parametrize('infile, expectedStdout', modifyName_params)
def test_modifyName(infile, expectedStdout, datadir, capsys):
    infile = datadir.join(infile)
    modifyName(infile)
    validateOutput(capsys.readouterr(), datadir, expectedStdout)
