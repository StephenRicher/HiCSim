#!/usr/bin/env bash

config=/home/stephen/phd/modelling/pipeline/config/local-config.yaml
snakemake --use-conda -kp --notemp --configfile "${config}" "${@}"
