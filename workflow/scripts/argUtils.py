#!/usr/bin/env python3

import sys
import re
import json
import logging
import argparse
import pandas as pd
import numpy as np


def setDefaults(parser):
    """ Set logging config and return args with associated function """

    args = parser.parse_args()
    logFormat = '%(asctime)s - %(levelname)s - %(funcName)s - %(message)s'
    try:
        logging.basicConfig(level=args.verbose, format=logFormat)
        del args.verbose
    except AttributeError:
        logging.basicConfig(format=logFormat)
        pass
    try:
        function = args.function
        del args.function
    except AttributeError:
        parser.print_help()
        sys.exit()

    return args, function


def createMainParent(verbose=True, version=None):
    """ Create parser of verbose/version to be added to parser/subparsers. """
    parent = argparse.ArgumentParser(add_help=False)
    if version:
        parent.add_argument(
            '--version', action='version', version=f'%(prog)s {version}')
    if verbose:
        parent.add_argument(
            '--verbose', action='store_const', const=logging.DEBUG,
            default=logging.INFO, help='verbose logging for debugging')
    return parent


def positiveInt(value):
    if (re.match(r'^[-+]?[0-9]+$', value)
            and int(value) >= 0):
        return int(value)
    else:
        raise argparse.ArgumentTypeError(
            f'{value} is not a valid positive integer.')
