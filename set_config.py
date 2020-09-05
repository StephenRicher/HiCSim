#!/usr/bin/env python3

import sys
import pandas as pd


class ConfigurationError(Exception):
    pass


def unset_defaults(default, outer_key=''):
    """ Return True if default or any value in default
        nested dictionary is an empty string.
    """
    out = 0
    if default == '':
        sys.stderr.write(
            f'\033[31mNo configuration provided for {outer_key} and '
             'no default available.\033[m\n')
        out += 1
    elif isinstance(default, dict):
        for key in default:
            full_key = f'{outer_key} : {key}' if outer_key else key
            out += unset_defaults(default[key], outer_key=full_key)
    return out


def set_default(config, default, key, outer_key):

    RC = 0
    if unset_defaults(default[key], outer_key=outer_key):
        RC = 1
    else:
        config[key] = default[key]
        sys.stderr.write(
            f'\033[33mNo configuration provided for {outer_key}'
            f' - setting to default: {config[key]}.\033[m\n')

    return RC


def set_config(config, default, outer_key=''):

    RC = 0
    for key in default:
        full_key = f'{outer_key} : {key}' if outer_key else key
        try:
            if isinstance(default[key], dict):
                set_config(config[key], default[key], outer_key=full_key)
            elif config[key] is None:
                RC += set_default(config, default, key, outer_key=full_key)
            else:
                sys.stderr.write(
                    f'\033[32mSetting {full_key} to: {config[key]}\033[m\n')
        except (KeyError, TypeError):
            RC += set_default(config, default, key, outer_key=full_key)

    if RC > 0:
        raise ConfigurationError(
            '\033[31mInvalid configuration setting.\033[m\n')

    return config


def read_paths(filepath):
    data = pd.read_table(filepath, names=['path'], squeeze=True)
    # Set index to start from 1 instead of 0
    data.index = [str(i+1) for i in range(len(data))]
    return data


def nBeads(file):
    """ Return number of beads in sequence file """
    nBeads = 0
    with open(file) as fh:
        for line in fh:
            nBeads += len(line.strip())
    return nBeads
