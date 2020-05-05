#!/usr/bin/env python3

import random
import tempfile
from set_config import set_config, read_paths

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

if not config:
    configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
# Note does not currently correctly deal with special case config['masking'] in set_config
default_config = {
    'name' :          'clustered_copolymer',
    'cluster':        False,
    'workdir':        workflow.basedir,
    'genome':         '',
    'build':          'genome',
    'ctcf':           None,
    'masking':
        [{'file':     ''         ,
          'character': ''         ,}],
    'region':         '',
    'chr':            '',
    'start':          '',
    'end':            '',
    'min_rep':        1,
    'bases_per_bead': 1000,
    'method':         'mean',
    'dpi':            600,
    'cmap':           'Reds',
    'transform':      'log10',
    'vmin':           0,
    'vmax':           None,
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
    'warm_up':        20000,
    'sim_time':       2000000,
    'delay':          10,
    'loop':           0,
    'tmpdir':         tempfile.gettempdir(),
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

rule all:
    input:
        [f'vmd/{NAME}.gif', f'qc/{NAME}-summed.png',
         f'sequence/{BUILD}-{REGION}.dat',
         f'genome/masked/{BUILD}-{REGION}-masked_all.fasta',
         expand('matrices/{name}-{rep}.h5', name = NAME, rep = REPS)]


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

if config['ctcf'] is not None:


    rule getOrientationCTCF:
        input:
            genome = rules.bgzipGenome.output,
            index = rules.indexGenome.output,
            ctcf = config['ctcf']
        output:
            'tracks/CTCF-{rep}.bed'
        log:
            'logs/getOrientationCTCF/{rep}.log'
        conda:
            f'{ENVS}/extract_ctcf_direction.yaml'
        shell:
            '{SCRIPTS}/extract_direction.py {input.genome} {input.ctcf} '
            '> {output} 2> {log}'


    rule sortBed:
        input:
            expand('tracks/CTCF-{rep}.bed',
                rep = range(1, len(config['ctcf']) + 1))
        output:
            'tracks/CTCF.sort.bed'
        group:
            'bedtools'
        log:
            'logs/sortBed.log'
        conda:
            f'{ENVS}/bedtools.yaml'
        shell:
            'bedtools sort -i <(cat {input}) > {output} 2> {log}'


    rule mergeBed:
        input:
            rules.sortBed.output
        output:
            'tracks/CTCF.merged.bed'
        group:
            'bedtools'
        log:
            'logs/mergeBed.log'
        conda:
            f'{ENVS}/bedtools.yaml'
        shell:
            'bedtools merge -s -c 4,5,6 -o count,median,distinct '
            '-i {input} > {output} 2> {log}'


    rule scaleBed:
        input:
            rules.mergeBed.output
        output:
            'tracks/CTCF.scaled.bed'
        group:
            'bedtools'
        log:
            'logs/scaleBed.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/scaleBedScore.py --log {input} > {output} 2> {log}'


    rule filterBedScore:
        input:
            rules.scaleBed.output
        output:
            'tracks/CTCF-sampled-{rep}.bed'
        params:
            rep = REPS,
            seed = lambda wildcards: SEEDS[REPS.index(int(wildcards.rep))]
        log:
            'logs/sampleBed/{rep}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/filterBedScore.py --seed {params.seed} '
            '{input} > {output} 2> {log}'


    rule splitOrientation:
        input:
            rules.filterBedScore.output
        output:
            forward = 'tracks/split/CTCF-forward-{rep}.bed',
            reverse = 'tracks/split/CTCF-reverse-{rep}.bed'
        params:
            min_rep = config['min_rep']
        group:
            'mask'
        log:
            'logs/splitOrientation/{rep}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/split_bed.py --min_rep {params.min_rep} '
            '--forward {output.forward} --reverse {output.reverse} '
            '{input} &> {log}'


def maskFastaInput(wildcards):
    if config['ctcf']:
        return [f'genome/{BUILD}.fa.gz'] + rules.splitOrientation.output
    else:
        return [f'genome/{BUILD}.fa.gz']


def getMasking(wc):
    """ Build masking command for maskFasta for track data """
    command = ''
    for entry in config['masking']:
        command += f'--bed {entry["file"]},{entry["character"]} '
    if config['ctcf']:
        command += (f'--bed tracks/filtered/split/CTCF-forward-{wc.rep}.bed,F '
                    f'--bed tracks/filtered/split/CTCF-reverse-{wc.rep}.bed,R ')
    return command


rule maskFasta:
    input:
        maskFastaInput
    output:
        temp(f'genome/masked/{BUILD}-{{rep}}-masked.fasta')
    params:
        tmpdir = config['tmpdir'],
        masking = getMasking
    group:
        'mask'
    log:
        'logs/maskFasta/{rep}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        '{SCRIPTS}/maskFasta.py --tmp {params.tmpdir} '
        '--genome <(zcat {input[0]}) {params.masking} '
        '> {output} 2> {log}'


rule indexMasked:
    input:
        rules.maskFasta.output
    output:
        f'{rules.maskFasta.output}.fai'
    group:
        'mask'
    log:
        'logs/indexMasked/{rep}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule extractRegion:
    input:
        fasta = rules.maskFasta.output,
        index = rules.indexMasked.output
    output:
        f'genome/masked/{BUILD}-{REGION}-{{rep}}-masked.fasta'
    group:
        'mask'
    params:
        region = f'{CHR}:{START}-{END}'
    log:
        f'logs/extractRegion/{BUILD}-{REGION}-{{rep}}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.fasta} {params.region} > {output} 2> {log}'


rule FastaToBeads:
    input:
        rules.extractRegion.output
    output:
        f'sequence/{BUILD}-{REGION}-{{rep}}-bead_sequence.txt'
    params:
        bases_per_bead = config['bases_per_bead']
    log:
        'logs/sequenceToBeads/{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/FastaToBeads.py --nbases {params.bases_per_bead} '
        '{input} > {output} 2> {log}'


rule BeadsToLammps:
    input:
        rules.FastaToBeads.output
    output:
        dat = f'sequence/{BUILD}-{REGION}-{{rep}}.dat',
        coeffs = f'sequence/{BUILD}-{REGION}-{{rep}}-ctcf_coeffs.txt',
    params:
        xlo = config['xlo'],
        xhi = config['xhi'],
        ylo = config['ylo'],
        yhi = config['yhi'],
        zlo = config['zlo'],
        zhi = config['zhi']
    log:
        f'logs/{BUILD}-{REGION}-{{rep}}-sequence_to_lammps.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_polymer.py model '
        '--xlo {params.xlo} --xhi {params.xhi} '
        '--ylo {params.ylo} --yhi {params.yhi} '
        '--zlo {params.zlo} --zhi {params.zhi} '
        '--ctcf --ctcfOut {output.coeffs} '
        '--ctcfCoeff 1.9,1.0,1.8 '
        '{input} > {output.dat} 2> {log}'


rule addCTCF:
    input:
        script = f'{SCRIPTS}/run-ctcf.lam',
        coeffs = rules.BeadsToLammps.output.coeffs
    output:
        'lammps/script/run-ctcf-{rep}.lam'
    log:
        'logs/addCTCF/{rep}.log'
    shell:
        "sed '/^#CTCF_COEFF/ r {input.coeffs}' {input.script} "
        "> {output} 2> {log}"


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


lmp_cmd = ('-var infile {input.data} '
           '-var outdir {params.outdir} '
           '-var name {NAME}-{wildcards.rep} '
           '-var warm_up {params.warm_up} '
           '-var sim_time {params.sim_time} '
           '-var timestep {params.timestep} '
           '-var seed {params.seed} '
           '-in {input.script} '
           '-log /dev/null &> {log}')
if config['cluster']:
    lmp_cmd = 'mpirun -n {threads} -x OMP_NUM_THREADS=1 lmp_mpi ' + lmp_cmd
else:
    lmp_cmd = 'lmp_serial ' + lmp_cmd


rule lammps:
    input:
        data = rules.BeadsToLammps.output.dat,
        script = rules.addCTCF.output
    output:
        warm_up = f'lammps/XYZ_{NAME}-{{rep}}-warm_up.xyz',
        proper_run = f'lammps/XYZ_{NAME}-{{rep}}-proper_run.xyz'
    params:
        outdir = directory(f'lammps'),
        timestep = config['timestep'],
        warm_up = config['warm_up'],
        sim_time = config['sim_time'],
        seed = lambda wildcards: SEEDS[REPS.index(int(wildcards.rep))]
    group:
        'lammps'
    threads:
        config['threads'] if config['cluster'] else 12
    log:
        'logs/lammps-{rep}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    shell:
        lmp_cmd


rule create_contact_matrix:
    input:
        rules.lammps.output.proper_run
    output:
        f'matrices/{NAME}-{{rep}}.npz'
    group:
        'lammps'
    log:
        'logs/contact_frequency-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/create_contact_matrix.py --outdata {output} {input} &> {log}'


rule matrix2homer:
    input:
        rules.create_contact_matrix.output
    output:
        f'matrices/{NAME}-{{rep}}.homer'
    params:
        chr = CHR,
        start = START,
        binsize = config['bases_per_bead']
    log:
        'logs/matrix2homer-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/npz2homer.py --chromosome {params.chr} '
        '--start {params.start} --binsize {params.binsize} '
        '{input} > {output} 2> {log}'


rule homerToH5:
    input:
        rules.matrix2homer.output
    output:
        f'matrices/{NAME}-{{rep}}.h5'
    log:
        'logs/homerToH5/{rep}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    threads:
        12
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 &> {log}'


rule average_matrices:
    input:
        expand('matrices/{name}-{rep}.npz',
            name=NAME, rep=REPS)
    output:
        f'qc/{NAME}-summed.npz'
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
        region = f'{START}-{END}',
        vmin = config['vmin'],
        vmax = config['vmax']
    log:
        'logs/plot_heatmap.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plot_heatmap.py --heatmap {output} --dpi {params.dpi} '
        '--vmin {params.vmin} --vmax {params.vmax} --region {params.region} '
        '--transform {params.transform} --cmap {params.cmap} {input} &> {log}'


rule mean_xyz:
    input:
        expand('lammps/XYZ_{name}-{rep}.xyz', name=NAME, rep=REPS)
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
        '{input} > {output} 2> {log}'


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
        '< {SCRIPTS}/vmd.tcl &> {log}'


def aggregateVMD(wildcards):
    """ Return all RGB files generated by vmd checkpoint. """

    dir = checkpoints.vmd.get(**wildcards).output[0]
    images = expand('{dir}/{i}.tga',
           i=glob_wildcards(os.path.join(dir, '{i}.tga')).i,
           dir=dir)
    return sorted(images)


rule create_gif:
    input:
        rules.vmd.output,
        images = aggregateVMD
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
        '{input.images} {output} &> {log}'
