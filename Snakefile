#!/usr/bin/env python3

import random
from set_config import set_config, read_paths

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
    'genome':       '',
    'build':        'genome',
    'ctcf':         '',
    'region':       '',
    'chr':          '',
    'start':        '',
    'end':          '',
    'min_rep':      1,
    'bases_per_bead': 1000,
    'method':       'median',
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

GENOME = config['genome']
BUILD = config['build']

REGION = config['region']
CHR = config['chr']
START = config['start']
END = config['end']
CTCF_DATA = read_paths(config['ctcf'])

rule all:
    input:
        [f'vmd/{NAME}.gif', f'qc/{NAME}-summed.png',
         f'sequence/{BUILD}-{REGION}.dat']
        #[f'qc/{NAME}-summed.png']


rule bgzip_genome:
    input:
        GENOME
    output:
        f'genome/{BUILD}.fa.gz'
    log:
        f'logs/bgzip_genome/{BUILD}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip > {output}) 2> {log}'


rule index_genome:
    input:
        rules.bgzip_genome.output
    output:
        multiext(f'{rules.bgzip_genome.output}','.fai', '.gzi')
    log:
        'logs/index_genome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule get_ctcf_orientation:
    input:
        genome = rules.bgzip_genome.output,
        index = rules.index_genome.output,
        ctcf = lambda wc: CTCF_DATA[wc.rep]
    output:
        'tracks/MCF-7_CTCF_{rep}.bed'
    params:
        chr = config['chr'],
    log:
        'logs/get_ctcf_orientation/{rep}.log'
    conda:
        f'{ENVS}/extract_ctcf_direction.yaml'
    shell:
        '{SCRIPTS}/extract_direction.py {input.genome} {input.ctcf} '
        '{params.chr} > {output} 2> {log}'


rule concatenate_bed:
    input:
        expand('tracks/MCF-7_CTCF_{rep}.bed', rep = CTCF_DATA.index)
    output:
        'tracks/MCF7_CTCF.bed'
    log:
        'logs/merge_bed.log'
    shell:
        'cat {input} > {output} 2> {log}'


rule bedtools_sort:
    input:
        rules.concatenate_bed.output
    output:
        'tracks/MCF7_CTCF.sort.bed'
    log:
        'logs/bedtools_sort.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools sort -i {input} > {output} 2> {log}'


rule bedtools_merge:
    input:
        rules.bedtools_sort.output
    output:
        'tracks/MCF7_CTCF.merged.bed'
    log:
        'logs/bedtools_merge.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools merge -s -c 4,5,6 -o count,median,distinct -i {input} '
        '> {output} 2> {log}'


rule subset_genome:
    input:
        genome = rules.bgzip_genome.output,
        indexes = rules.index_genome.output
    output:
        f'genome/{BUILD}-{CHR}.fa'
    params:
        chr = config['chr'],
    log:
        'logs/subset_genome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.genome} {params.chr} '
        '| sed "1s/://" > {output} 2> {log}'


rule index_subset_genome:
    input:
        rules.subset_genome.output
    output:
        f'{rules.subset_genome.output}.fai'
    log:
        'logs/index_subset_genome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule filter_and_split_orientation:
    input:
        rules.bedtools_merge.output
    output:
        forward = 'tracks/MCF7_CTCF-forward.bed',
        reverse = 'tracks/MCF7_CTCF-reverse.bed'
    params:
        min_rep = config['min_rep']
    log:
        'logs/split_bed.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/split_bed.py --min_rep {params.min_rep} '
        '--forward {output.forward} --reverse {output.reverse} {input} &> {log}'


rule merge_filtered:
    input:
        rules.filter_and_split_orientation.output
    output:
        'tracks/MCF7_CTCF-filtered.bed'
    log:
        'logs/merge_filtered.log'
    shell:
        'cat {input} > {output}'


rule sort_filtered:
    input:
        rules.merge_filtered.output
    output:
        'tracks/MCF7_CTCF-filtered.sort.bed'
    log:
        'logs/sort_filtered.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools sort -i {input} > {output} 2> {log}'


rule get_chrom_sizes:
    input:
        rules.index_subset_genome.output
    output:
        f'genome/chrom_sizes/{BUILD}-{CHR}.size'
    params:
        chr = config['chr']
    log:
        f'logs/get_chrom_sizes/{BUILD}-{CHR}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


rule bedtools_complement:
    input:
        sizes = rules.get_chrom_sizes.output,
        bed = rules.sort_filtered.output
    output:
        f'genome/{BUILD}-{CHR}-complement.bed'
    log:
        f'logs/bedtools_complement/{BUILD}-{CHR}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools complement -i {input.bed} -g {input.sizes} '
        '> {output} 2> {log}'


rule mask_forward_ctcf:
    input:
        bed = rules.filter_and_split_orientation.output.forward,
        genome = rules.subset_genome.output
    output:
        f'genome/masked/{BUILD}-{CHR}-masked_F.fasta'
    params:
        mc = 'F'
    log:
        'logs/mask_forward_ctcf.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools maskfasta -bed {input.bed} -mc {params.mc} '
        '-fi {input.genome} -fo {output} &> {log}'


rule mask_reverse_ctcf:
    input:
        bed = rules.filter_and_split_orientation.output.reverse,
        masked = rules.mask_forward_ctcf.output
    output:
        f'genome/masked/{BUILD}-{CHR}-masked_RF.fasta'
    params:
        mc = 'R'
    log:
        'logs/mask_reverse_ctcf.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools maskfasta -bed {input.bed} -mc {params.mc} '
        '-fi {input.masked} -fo {output} &> {log}'


rule masked_other:
    input:
        bed = rules.bedtools_complement.output,
        masked = rules.mask_reverse_ctcf.output
    output:
        f'genome/masked/{BUILD}-{CHR}-masked_all.fasta'
    params:
        mc = 'N'
    log:
        'logs/mask_other.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools maskfasta -bed {input.bed} -mc {params.mc} '
        '-fi {input.masked} -fo {output} &> {log}'

rule index_masked:
    input:
        rules.masked_other.output
    output:
        f'{rules.masked_other.output}.fai'
    log:
        'logs/index_masked.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule subset_masked:
    input:
        fasta = rules.masked_other.output,
        index = rules.index_masked.output
    output:
        f'genome/masked/{BUILD}-{REGION}-masked_all.fasta'
    params:
        region = f'{CHR}:{START}-{END}'
    log:
        f'logs/{BUILD}-{REGION}-subset_masked.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.fasta} {params.region} > {output} 2> {log}'


rule compress_sequence:
    input:
        rules.subset_masked.output
    output:
        f'sequence/{BUILD}-{REGION}-bead_sequence.txt'
    params:
        bases_per_bead = config['bases_per_bead']
    log:
        'logs/compress_sequence.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/compress_sequence.py --nbases {params.bases_per_bead} '
        '{input} > {output} 2> {log}'


rule sequence_to_lammps:
    input:
        rules.compress_sequence.output
    output:
        f'sequence/{BUILD}-{REGION}.dat'
    params:
        xlo = config['xlo'],
        xhi = config['xhi'],
        ylo = config['ylo'],
        yhi = config['yhi'],
        zlo = config['zlo'],
        zhi = config['zhi']
    log:
        f'logs/{BUILD}-{REGION}-sequence_to_lammps.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_polymer.py model '
        '--xlo {params.xlo} --xhi {params.xhi} '
        '--ylo {params.ylo} --yhi {params.yhi} '
        '--zlo {params.zlo} --zhi {params.zhi} '
        '{input} > {output} 2> {log}'

rule generate_random_polymer:
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
        rules.sequence_to_lammps.output
    output:
        warm_up = f'lammps/XYZ_{NAME}-{{rep}}-warm_up.xyz',
        proper_run = f'lammps/XYZ_{NAME}-{{rep}}-proper_run.xyz'
    params:
        rep = REPS,
        outdir = directory(f'lammps'),
        timestep = config['timestep'],
        seed = lambda wildcards: SEEDS[REPS.index(int(wildcards.rep))]
    threads:
        config['threads']
    log:
        'logs/lammps-{rep}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    shell:
        'mpirun -n {threads} -x OMP_NUM_THREADS={threads} '
        'lmp_mpi -var infile {input} '
        '-var outdir {params.outdir} '
        '-var name {NAME}-{wildcards.rep} '
        '-var timestep {params.timestep} '
        '-var seed {params.seed} '
        '-in {SCRIPTS}/run-ctcf.lam '
        '-log /dev/null &> {log}'


rule create_contact_matrix:
    input:
        rules.lammps.output.proper_run
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


rule average_matrices:
    input:
        expand('qc/{name}-{rep}.txt.gz',
            name=NAME, rep=REPS)
    output:
        f'qc/{NAME}-summed.txt.gz'
    params:
        method = config['method']
    log:
        'logs/sum_matrices.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/average_matrices.py --out {output} '
        '--method {params.method} {input} &> {log}'


rule plot_heatmap:
    input:
        rules.average_matrices.output
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
        f'lammps/XYZ_{NAME}-1-proper_run.xyz'
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
        'vmd -eofexit -nt -dispdev text -args {input} {output} '
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
