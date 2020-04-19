#!/usr/bin/env python3

import random
from set_config import set_config, read_paths

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

if not config:
    configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'name' :          'clustered_copolymer',
    'workdir':        workflow.basedir,
    'genome':         '',
    'build':          'genome',
    'ctcf':           '',
    'region':         '',
    'chr':            '',
    'start':          '',
    'end':            '',
    'min_rep':        1,
    'bases_per_bead': 1000,
    'method':         'mean',
    'dpi':            600,
    'cmap':           'YlGn',
    'transform':      'obsexp',
    'vmin':           0,
    'vmax':           2,
    'n_molecules':    1000,
    'reps':           5,
    'threads':        1,
    'seed':           42,
    'n_types':        4,
    'n_clusters':     20,
    'xlo':            -100,
    'xhi':            100,
    'ylo':            -100,
    'yhi':            100,
    'zlo':            -100,
    'zhi':            100,
    'window_size':    10,
    'overlap':        2,
    'timestep':       1000,
    'delay':          10,
    'loop':           0
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


rule bgzipGenome:
    input:
        GENOME
    output:
        f'genome/{BUILD}.fa.gz'
    log:
        f'logs/bgzipGenome/{BUILD}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        '(zcat -f {input} | bgzip > {output}) 2> {log}'


rule indexGenome:
    input:
        rules.bgzipGenome.output
    output:
        multiext(f'{rules.bgzipGenome.output}', '.fai', '.gzi')
    log:
        'logs/indexGenome.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule get_ctcf_orientation:
    input:
        genome = rules.bgzipGenome.output,
        index = rules.indexGenome.output,
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


rule catReplicates:
    input:
        expand('tracks/MCF-7_CTCF_{rep}.bed', rep = CTCF_DATA.index)
    output:
        'tracks/MCF7_CTCF.bed'
    group:
        'bedtools'
    log:
        'logs/merge_bed.log'
    shell:
        'cat {input} > {output} 2> {log}'


rule sortBed:
    input:
        rules.catReplicates.output
    output:
        'tracks/MCF7_CTCF.sort.bed'
    group:
        'bedtools'
    log:
        'logs/sortBed.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools sort -i {input} > {output} 2> {log}'


rule mergeBed:
    input:
        rules.sortBed.output
    output:
        'tracks/MCF7_CTCF.merged.bed'
    group:
        'bedtools'
    log:
        'logs/mergeBed.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools merge -s -c 4,5,6 -o count,median,distinct -i {input} '
        '> {output} 2> {log}'


rule splitOrientationFilterBed:
    input:
        rules.mergeBed.output
    output:
        forward = 'tracks/MCF7_CTCF-forward.bed',
        reverse = 'tracks/MCF7_CTCF-reverse.bed'
    params:
        min_rep = config['min_rep']
    group:
        'bedtools'
    log:
        'logs/split_bed.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/split_bed.py --min_rep {params.min_rep} '
        '--forward {output.forward} --reverse {output.reverse} {input} &> {log}'


rule catBed:
    input:
        rules.splitOrientationFilterBed.output
    output:
        'tracks/MCF7_CTCF-filtered.bed'
    group:
        'bedtools'
    log:
        'logs/catBed.log'
    shell:
        'cat {input} > {output}'


rule sortFilteredBed:
    input:
        rules.catBed.output
    output:
        'tracks/MCF7_CTCF-filtered.sort.bed'
    group:
        'bedtools'
    log:
        'logs/sortFilteredBed.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools sort -i {input} > {output} 2> {log}'


rule extractChrom:
    input:
        genome = rules.bgzipGenome.output,
        indexes = rules.indexGenome.output
    output:
        f'genome/{BUILD}-{CHR}.fa'
    params:
        chr = config['chr'],
    group:
        'extractChrom'
    log:
        'logs/extractChrom.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.genome} {params.chr} '
        '| sed "1s/://" > {output} 2> {log}'


rule indexChrom:
    input:
        rules.extractChrom.output
    output:
        f'{rules.extractChrom.output}.fai'
    group:
        'extractChrom'
    log:
        'logs/indexChrom.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule getChromSize:
    input:
        rules.indexChrom.output
    output:
        f'genome/chrom_sizes/{BUILD}-{CHR}.size'
    params:
        chr = config['chr']
    group:
        'extractChrom'
    log:
        f'logs/getChromSize/{BUILD}-{CHR}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


rule complementBed:
    input:
        sizes = rules.getChromSize.output,
        bed = rules.sortFilteredBed.output
    output:
        f'genome/{BUILD}-{CHR}-complement.bed'
    log:
        f'logs/complementBed/{BUILD}-{CHR}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools complement -i {input.bed} -g {input.sizes} '
        '> {output} 2> {log}'


rule maskForwardCTCF:
    input:
        bed = rules.splitOrientationFilterBed.output.forward,
        genome = rules.extractChrom.output
    output:
        f'genome/masked/{BUILD}-{CHR}-masked_F.fasta'
    params:
        mc = 'F'
    log:
        'logs/maskForwardCTCF.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools maskfasta -bed {input.bed} -mc {params.mc} '
        '-fi {input.genome} -fo {output} &> {log}'


rule maskReverseCTCF:
    input:
        bed = rules.splitOrientationFilterBed.output.reverse,
        masked = rules.maskForwardCTCF.output
    output:
        f'genome/masked/{BUILD}-{CHR}-masked_RF.fasta'
    params:
        mc = 'R'
    group:
        'mask'
    log:
        'logs/maskReverseCTCF.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools maskfasta -bed {input.bed} -mc {params.mc} '
        '-fi {input.masked} -fo {output} &> {log}'


rule maskOther:
    input:
        bed = rules.complementBed.output,
        masked = rules.maskReverseCTCF.output
    output:
        f'genome/masked/{BUILD}-{CHR}-masked_all.fasta'
    params:
        mc = 'N'
    group:
        'mask'
    log:
        'logs/mask_other.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        'bedtools maskfasta -bed {input.bed} -mc {params.mc} '
        '-fi {input.masked} -fo {output} &> {log}'

rule indexMasked:
    input:
        rules.maskOther.output
    output:
        f'{rules.maskOther.output}.fai'
    group:
        'mask'
    log:
        'logs/indexMasked.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule extractRegion:
    input:
        fasta = rules.maskOther.output,
        index = rules.indexMasked.output
    output:
        f'genome/masked/{BUILD}-{REGION}-masked_all.fasta'
    group:
        'mask'
    params:
        region = f'{CHR}:{START}-{END}'
    log:
        f'logs/{BUILD}-{REGION}-extractRegion.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.fasta} {params.region} > {output} 2> {log}'


rule sequenceToBeads:
    input:
        rules.extractRegion.output
    output:
        f'sequence/{BUILD}-{REGION}-bead_sequence.txt'
    params:
        bases_per_bead = config['bases_per_bead']
    log:
        'logs/sequenceToBeads.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/sequenceToBeads.py --nbases {params.bases_per_bead} '
        '{input} > {output} 2> {log}'


rule sequence_to_lammps:
    input:
        rules.sequenceToBeads.output
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
    group:
        'lammps'
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
    group:
        'lammps'
    log:
        'logs/contact_frequency-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/create_contact_matrix.py --outdata {output} {input} &> {log}'


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
    params:
        dpi = config['dpi'],
        cmap = config['cmap'],
        transform = config['transform'],
        vmin = config['vmin'],
        vmax = config['vmax']
    log:
        'logs/plot_heatmap.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plot_heatmap.py --heatmap {output} --dpi {params.dpi} '
        '--vmin {params.vmin} --vmax {params.vmax} '
        '--transform {params.transform} --cmap {params.cmap} {input} &> {log}'


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
