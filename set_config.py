#!/usr/bin/env python3

import sys
import pandas as pd


class ConfigurationError(Exception):
    pass


def set_config(config, default_config):

    RC = 0

    for key in default_config:
        try:
            config[key]
            sys.stderr.write(f'\033[32mSetting {key} to: {config[key]}\n')
        except KeyError:
            if default_config[key] == '':
                sys.stderr.write(
                    f'\033[31mNo configuration provided for {key} and '
                     'no default available.\n')
                RC = 1
            else:
                config[key] = default_config[key]
                sys.stderr.write(
                    f'\033[33mNo configuration provided for {key}.\n')
                sys.stderr.write(
                    f'\033[33mSetting {key} to default: {config[key]}.\n')

    if RC == 1:
        raise ConfigurationError(
            '\033[31mInvalid configuration setting.\033[m\n')

    sys.stderr.write('\033[m')
    return config


def read_paths(filepath):
    data = pd.read_table(filepath, names=['path'], squeeze=True)
    # Set index to start from 1 instead of 0
    data.index = [str(i+1) for i in range(len(data))]
    return data
