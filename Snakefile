#!/usr/bin/env python3

import os
import random
import tempfile
from set_config import set_config, read_paths

wildcard_constraints:
    all = r'[^\/]+',

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
    'cluster':        False,
    'workdir':        workflow.basedir,
    'genome':         '',
    'build':          'genome',
    'ctcf':           None,
    'masking':        {},
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
    'HiC':            {'matrix' : None,
                       'binsize': None,
                       'log' :    True},
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
    'coeffs':         ''
}
config = set_config(config, default_config)

workdir : config['workdir']
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

## Replace file name to generate dictionary of modified BEDS
squeezed_bed = {}
for file, character in config['masking'].items():
    if file in squeezed_bed:
        sys.exit(f'Duplicate file {file} in config["masking"].')
    squeezed_file = f'{os.path.basename(file)}'
    squeezed_bed[squeezed_file] = {'source' : file, 'character' : character}

# ADD TO CONFIG
NMONOMERS = 200
# Also see pairCoeffs in Beads2Lammps

BINSIZE = int(config['HiC']['binsize'])
if BINSIZE:
    MERGEBINS, REMAINDER = divmod(BINSIZE, config['bases_per_bead'])
    if REMAINDER:
        sys.exit(f'Binsize {config["HiC"]["binsize"]} is not divisible by '
            f'bases per bead {config["bases_per_bead"]}')
else:
    MERGEBINS = 1


rule all:
    input:
        [f'vmd/simulation.gif', 'merged/simulation.png']


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


rule getChromSizes:
    input:
        rules.indexGenome.output
    output:
        f'genome/chrom_sizes/{BUILD}.chrom.sizes'
    log:
        f'logs/getChromSizes/{BUILD}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


rule extractChrom:
    input:
        fasta = rules.bgzipGenome.output,
        index = rules.indexGenome.output
    output:
        f'genome/{BUILD}-{CHR}.fa'
    params:
        chrom = CHR
    log:
        'logs/extractChrom.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.fasta} {params.chrom} > {output} 2> {log}'


rule indexChrom:
    input:
        rules.extractChrom.output
    output:
        multiext(f'{rules.extractChrom.output}', '.fai')
    log:
        'logs/indexChrom.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


if config['ctcf'] is not None:

    rule getOrientationCTCF:
        input:
            genome = rules.extractChrom.output,
            index = rules.indexChrom.output,
            ctcf = lambda wc: config['ctcf'][int(wc.rep)]
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
                rep = range(0, len(config['ctcf'])))
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


    rule squeezeCTCF:
        input:
            rules.mergeBed.output
        output:
            'tracks/CTCF-squeezed.bed'
        params:
            length = 20
        log:
            'logs/squeezeCTCF.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/squeezeBed.py {params.length} {input} '
            '> {output} 2> {log}'


    rule scaleBed:
        input:
            rules.squeezeCTCF.output
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
            'replicates/{rep}/tracks/CTCF-sampled-{rep}.bed'
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
            forward = 'replicates/{rep}/tracks/CTCF-forward-{rep}.bed',
            reverse = 'replicates/{rep}/tracks/CTCF-reverse-{rep}.bed'
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


rule squeezeBed:
    input:
        lambda wc: squeezed_bed[wc.squeezed_bed]['source']
    output:
        'tracks/squeezed/{squeezed_bed}'
    params:
        length = 20
    log:
        'logs/squeezeBed/{squeezed_bed}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/squeezeBed.py {params.length} {input} '
        '> {output} 2> {log}'


def getMasking(wc):
    """ Build masking command for maskFasta for track data """
    command = ''
    for bed in squeezed_bed:
        character = squeezed_bed[bed]['character']
        command += f'--bed tracks/squeezed/{bed},{character} '
    if config['ctcf']:
        command += (f'--bed replicates/{wc.rep}/tracks/CTCF-forward-{wc.rep}.bed,F '
                    f'--bed replicates/{wc.rep}/tracks/CTCF-reverse-{wc.rep}.bed,R ')
    return command


rule maskFasta:
    input:
        expand('tracks/squeezed/{squeezed_bed}',
            squeezed_bed = squeezed_bed.keys()),
        rules.splitOrientation.output,
        genome = rules.extractChrom.output
    output:
        pipe(f'replicates/{{rep}}/genome/{BUILD}-masked-{{rep}}.fa')
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
        '--genome {input.genome} {params.masking} '
        '> {output} 2> {log}'

rule bgzipMasked:
    input:
        rules.maskFasta.output
    output:
        f'{rules.maskFasta.output}.gz'
    group:
        'mask'
    log:
        f'logs/bgzipMasked/{BUILD}-{{rep}}.log'
    conda:
        f'{ENVS}/tabix.yaml'
    shell:
        'bgzip -c {input} > {output} 2> {log}'

rule indexMasked:
    input:
        rules.bgzipMasked.output
    output:
        multiext(f'{rules.bgzipMasked.output}', '.fai', '.gzi')
    group:
        'mask'
    log:
        'logs/indexMasked/{rep}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'


rule extractRegion:
    input:
        fasta = rules.bgzipMasked.output,
        index = rules.indexMasked.output
    output:
        pipe(f'replicates/{{rep}}/genome/{BUILD}-{REGION}-masked-{{rep}}.fa')
    group:
        'convert2Beads'
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
        f'replicates/{{rep}}/sequence/{BUILD}-{REGION}-beads-{{rep}}.txt'
    group:
        'convert2Beads'
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
        dat = f'replicates/{{rep}}/lammps/config/{BUILD}-{REGION}-{{rep}}.dat',
        coeffs = f'replicates/{{rep}}/lammps/config/{BUILD}-{REGION}-{{rep}}-coeffs.txt',
        groups = f'replicates/{{rep}}/lammps/config/{BUILD}-{REGION}-{{rep}}-groups.txt'
    params:
        nmonomers = NMONOMERS,
        coeffs = config['coeffs'],
        basesPerBead = config['bases_per_bead'],
        seed = lambda wc: SEEDS[REPS.index(int(wc.rep))],
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
        '{SCRIPTS}/generate_polymer.py model --seed {params.seed} '
        '--xlo {params.xlo} --xhi {params.xhi} '
        '--ylo {params.ylo} --yhi {params.yhi} '
        '--zlo {params.zlo} --zhi {params.zhi} '
        '--ctcf --coeffOut {output.coeffs} '
        '--groupOut {output.groups} '
        '--monomer {params.nmonomers},T '
        '--pairCoeffs {params.coeffs} '
        '--basesPerBead {params.basesPerBead} '
        '{input} > {output.dat} 2> {log}'


rule addCTCF:
    input:
        script = f'{SCRIPTS}/run-ctcf.lam',
        coeffs = rules.BeadsToLammps.output.coeffs,
        groups = rules.BeadsToLammps.output.groups
    output:
        'replicates/{rep}/lammps/config/run-ctcf-{rep}.lam'
    log:
        'logs/addCTCF/{rep}.log'
    shell:
        "sed '/^#CTCF_COEFF/ r {input.coeffs}' {input.script} "
        "| sed '/^#GROUPS/ r {input.groups}' > {output} 2> {log}"


lmp_cmd = ('-var infile {input.data} '
           '-var radius_gyration {output.radius_gyration} '
           '-var warm_up {output.warm_up} '
           '-var warm_up_time {params.warm_up_time} '
           '-var restart {params.restart} '
           '-var restart_time {params.restart_time} '
           '-var sim {output.simulation} '
           '-var sim_time {params.sim_time} '
           '-var complete {output.complete} '
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
        warm_up = 'replicates/{rep}/lammps/warm_up-{rep}.xyz.gz',
        simulation = 'replicates/{rep}/lammps/simulation-{rep}.xyz.gz',
        complete = 'replicates/{rep}/lammps/complete-{rep}.xyz.gz',
        radius_gyration = 'replicates/{rep}/lammps/radius_of_gyration-{rep}.txt',
        restart = directory('replicates/{rep}/lammps/restart/')
    params:
        timestep = config['timestep'],
        warm_up_time = config['warm_up'],
        sim_time = config['sim_time'],
        restart_time = int(config['timestep']) * 100, # Use 0 for no restart file
        restart = lambda wc: f'replicates/{wc.rep}/lammps/restart/Restart-{wc.rep}',
        seed = lambda wc: SEEDS[REPS.index(int(wc.rep))]
    group:
        'lammps'
    threads:
        config['threads'] if config['cluster'] else 1
    log:
        'logs/lammps-{rep}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    shell:
        lmp_cmd


rule create_contact_matrix:
    input:
        rules.lammps.output.simulation
    output:
        'replicates/{rep}/matrices/{rep}.npz'
    group:
        'lammps'
    log:
        'logs/contact_frequency-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/create_contact_matrix.py --outdata {output} '
        '<(zcat {input}) &> {log}'


rule mergeReplicates:
    input:
        expand('replicates/{rep}/matrices/{rep}.npz', rep=REPS)
    output:
        'merged/merged.npz'
    params:
        method = config['method']
    log:
        'logs/mergeReplicates.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/average_matrices.py --out {output} '
        '--method {params.method} {input} &> {log}'


rule matrix2homer:
    input:
        rules.mergeReplicates.output
    output:
        'merged/merged.homer'
    params:
        chr = CHR,
        start = START,
        binsize = config['bases_per_bead']
    log:
        'logs/matrix2homer.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/npz2homer.py --chromosome {params.chr} '
        '--start {params.start} --binsize {params.binsize} '
        '{input} > {output} 2> {log}'


rule homer2H5:
    input:
        rules.matrix2homer.output
    output:
        'merged/merged.h5'
    log:
        'logs/homer2H5.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    threads:
        12
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 &> {log}'


rule mergeBins:
    input:
        rules.homer2H5.output
    output:
        f'merged/merged-{BINSIZE}.h5'
    params:
        nbins = MERGEBINS
    log:
        'logs/mergeBins.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log}'


rule matrix2pre:
    input:
        rules.mergeReplicates.output
    output:
        'merged/merged.pre.tsv'
    params:
        chr = CHR,
        start = START,
        binsize = config['bases_per_bead']
    log:
        'logs/matrix2pre.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/npz2pre.py --chromosome {params.chr} '
        '--start {params.start} --binsize {params.binsize} '
        '{input} > {output} 2> {log}'


rule juicerPre:
    input:
        tsv = rules.matrix2pre.output,
        chrom_sizes = rules.getChromSizes.output
    output:
        'merged/merged.hic'
    log:
        'logs/juicerPre.log'
    params:
        chr = CHR,
        resolutions = '500,1000,1500'
    resources:
         mem_mb = 16000
    threads:
        12
    conda:
        f'{ENVS}/openjdk.yaml'
    shell:
        'java -Xmx{resources.mem_mb}m '
        '-jar {SCRIPTS}/juicer_tools_1.14.08.jar pre '
        '-c {params.chr} -r {params.resolutions} '
        '{input.tsv} {output} {input.chrom_sizes} &> {log}'


def getHiCconfig(wc):
    """ Build config command for real HiC data """
    command = ''
    if config['HiC']['matrix']:
        command += f'--flip --matrix2 {config["HiC"]["matrix"]} '
    if config['HiC']['log']:
        command += '--log_matrix2'
    return command


rule createConfig:
    input:
        matrix = rules.mergeBins.output,
        ctcf_orientation = rules.scaleBed.output,
        genes = config['genes']
    output:
        'merged/configs.ini'
    conda:
        f'{ENVS}/python3.yaml'
    params:
        depth = int((END - START) / 2),
        hicConfig = getHiCconfig
    log:
        'logs/createConfig.log'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.matrix} '
        '--ctcf_orientation {input.ctcf_orientation} --log '
        '--genes {input.genes} {params.hicConfig} '
        '--depth {params.depth} > {output} 2> {log}'


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        'merged/simulation.png'
    params:
        region = f'{CHR}:{START}-{END}',
        title = f'"{REGION} : {CHR}:{START}-{END}"',
        dpi = 600
    log:
        'logs/plotHiC.log'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    shell:
        'pyGenomeTracks --tracks {input} '
        '--region {params.region} '
        '--outFileName {output} '
        '--title {params.title} '
        '--dpi {params.dpi} &> {log}'


rule smooth_xyz:
    input:
        'replicates/1/lammps/complete-1.xyz.gz'
    output:
        'replicates/1/lammps/complete-smooth-1.xyz'
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
        '<(zcat {input}) > {output} 2> {log}'


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
        f'vmd/simulation.gif'
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
