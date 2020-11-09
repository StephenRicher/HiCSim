#!/usr/bin/env python

import os
import sys
import pytest
from testUtils import datadir, validateOutput
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from subsetATAC import *


main_params = (
    [('test1-100.in', {'chr': '1', 'start': 0, 'end' : 1000}, 100, 'test1-100.out')])
@pytest.mark.parametrize(
    'infile, region, nbases, expectedStdout', main_params)
def test_main(infile, region, nbases, expectedStdout, datadir, capsys):
    infile = datadir.join(infile)
    main(infile, region, nbases)
    validateOutput(capsys.readouterr(), datadir, expectedStdout)


validCoordinates_params = (
    [(0, 1000, 10, True ),
     (1, 1000, 10, False),
     (1,    1,  5, False)])
@pytest.mark.parametrize(
    'start, end, nbases, out', validCoordinates_params)
def test_validCoordinates(start, end, nbases, out):
    assert out == validCoordinates(start, end, nbases)
