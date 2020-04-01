#!/usr/bin/env python3

import random
from set_config import set_config

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'name' :        'clustered_copolymer',
    'workdir':      workflow.basedir,
    'n_molecules':  1000,
    'reps':         5,
    'threads':      1,
    'seed':         42,
    'n_types':      4,
    'n_clusters':   20,
    'xlo':          -100,
    'xhi':          100,
    'ylo':          -100,
    'yhi':          100,
    'zlo':          -100,
    'zhi':          100,
    'window_size':  10,
    'overlap':      2,
    'timestep':     1000,
    'delay':        10,
    'loop':         0
}
config = set_config(config, default_config)

workdir : config['workdir']
# Define global configuration variables
NAME = config['name']
# Define list of N reps from 1 to N
REPS = list(range(1, config['reps'] + 1))
# Define random seed for each rep
random.seed(config['seed'])
SEEDS = [random.randint(1, (2**16) - 1) for i in REPS]

rule all:
    input:
        [f'vmd/{NAME}.gif', f'qc/{NAME}-summed.png']
        #[f'qc/{NAME}-summed.png']


rule generate_polymer:
    output:
        f'polymer/poly.n{config["n_molecules"]}.dat'
    params:
        n_molecules = config['n_molecules'],
        n_clusters = config['n_clusters'],
        n_types = config['n_types'],
        xlo = config['xlo'],
        xhi = config['xhi'],
        ylo = config['ylo'],
        yhi = config['yhi'],
        zlo = config['zlo'],
        zhi = config['zhi'],
    log:
        'logs/generate_polymer.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_polymer.py '
            '--n_molecules {params.n_molecules} '
            '--n_clusters {params.n_clusters} '
            '--n_types {params.n_types} '
            '--xlo {params.xlo} --xhi {params.xhi} '
            '--ylo {params.ylo} --yhi {params.yhi} '
            '--zlo {params.zlo} --zhi {params.zhi} '
        '> {output} 2> {log}'


rule lammps:
    input:
        rules.generate_polymer.output
    output:
        f'lammps/XYZ_{NAME}-{{rep}}.xyz'
    params:
        rep = REPS,
        outdir = directory(f'lammps'),
        timestep = config['timestep'],
        seed = lambda wildcards: SEEDS[REPS.index(int(wildcards.rep))]
    log:
        'logs/lammps-{rep}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    shell:
        'lmp_serial '
            '-var infile {input} '
            '-var outdir {params.outdir} '
            '-var name {NAME}-{wildcards.rep} '
            '-var timestep {params.timestep} '
            '-var seed {params.seed} '
            '-in {SCRIPTS}/run.lam '
            '-log /dev/null '
        '&> {log}'


rule create_contact_matrix:
    input:
        rules.lammps.output
    output:
        f'qc/{NAME}-{{rep}}.txt.gz'
    params:
        xsize = config['xhi'] - config['xlo'],
        ysize = config['yhi'] - config['ylo'],
        zsize = config['zhi'] - config['zlo']
    threads:
        config['threads']
    log:
        'logs/contact_frequency-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/create_contact_matrix.py '
            '--threads {threads} '
            '--xsize {params.xsize} '
            '--ysize {params.ysize} '
            '--zsize {params.zsize} '
            '--outdata {output} {input} '
        '&> {log}'


rule sum_matrices:
    input:
        expand('qc/{name}-{rep}.txt.gz',
            name=NAME, rep=REPS)
    output:
        f'qc/{NAME}-summed.txt.gz'
    log:
        'logs/sum_matrices.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/sum_matrices.py '
            '--out {output} {input} '
        '&> {log}'


rule plot_heatmap:
    input:
        rules.sum_matrices.output
    output:
        f'qc/{NAME}-summed.png'
    log:
        'logs/plot_heatmap.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plot_heatmap.py '
            '--heatmap {output} {input} '
        '&> {log}'


rule mean_xyz:
    input:
        expand('lammps/XYZ_{name}-{rep}.xyz',
            name=NAME, rep=REPS)
    output:
        f'lammps/XYZ_{NAME}-mean.xyz'
    log:
        'logs/mean_xyz.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mean_xyz.py --verbose {input} > {output} 2> {log}'


rule smooth_xyz:
    input:
        rules.mean_xyz.output
    output:
        f'lammps/XYZ_{NAME}-smooth.xyz'
    params:
        window_size = config['window_size'],
        overlap = config['overlap']
    log:
        'logs/smooth_xyz.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/smooth_xyz.py '
            '--size {params.window_size} '
            '--overlap {params.overlap} '
            '{input} '
        '> {output} 2> {log}'


checkpoint vmd:
    input:
        rules.smooth_xyz.output
    output:
        directory(f'vmd/sequence')
    log:
        'logs/vmd.log'
    conda:
        f'{ENVS}/vmd.yaml'
    shell:
        'vmd -eofexit -nt -displaydev text -args {input} {output} '
            '< {SCRIPTS}/vmd.tcl '
        '&> {log}'


def aggregate_vmd(wildcards):
    """ Return all RGB files generated by vmd checkpoint. """

    dir = checkpoints.vmd.get(**wildcards).output[0]
    images = expand('{dir}/{i}.tga',
           i=glob_wildcards(os.path.join(dir, '{i}.tga')).i,
           dir=dir)
    return sorted(images)


rule create_gif:
    input:
        rules.vmd.output,
        images = aggregate_vmd
    output:
        f'vmd/{NAME}.gif'
    params:
        delay = config['delay'],
        loop = config['loop']
    log:
        'logs/create_gif.log'
    conda:
        f'{ENVS}/imagemagick.yaml'
    shell:
        'convert -delay {params.delay} -loop {params.loop} '
            '{input.images} {output} '
        '&> {log}'
