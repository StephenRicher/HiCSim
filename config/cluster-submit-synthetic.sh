#!/usr/bin/env bash

config=/home/u/sr467/scratch/projects/MolecularDynamics/config/config-synthetic.yaml
cluster_config=/home/u/sr467/scratch/projects/MolecularDynamics/config/cluster-config.yaml
conda_prefix=/home/u/sr467/scratch/projects/MolecularDynamics/.snakemake/

snakemake --profile slurm --use-conda -kp --notemp \
          --cluster-config "${cluster_config}" --configfile "${config}" \
          --conda-prefix "${conda_prefix}" "${@}"
