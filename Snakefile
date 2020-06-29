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
    'monomers':       100,
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

BINSIZE = int(config['HiC']['binsize'])
if BINSIZE:
    MERGEBINS, REMAINDER = divmod(BINSIZE, config['bases_per_bead'])
    if REMAINDER:
        sys.exit(
            f'Binsize {config["HiC"]["binsize"]} is not divisible by '
            f'bases per bead {config["bases_per_bead"]}')
else:
    MERGEBINS = 1
    BINSIZE = config['bases_per_bead']


track_data = {}
for file, character in config['masking'].items():
    track_file = f'{os.path.basename(file)}'
    track_data[track_file] = {'source' : file, 'character' : character}

rule all:
    input:
        [expand('{region}/{nbases}/vmd/simulation.gif',
            nbases=config['bases_per_bead'], region=REGION),
         expand('{region}/{nbases}/merged/simulation-{binsize}.png',
            nbases=config['bases_per_bead'], region=REGION, binsize=BINSIZE)]


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
        fai = f'{rules.bgzipGenome.output}.fai',
        gzi = f'{rules.bgzipGenome.output}.gzi'
    log:
        f'logs/indexGenome/{BUILD}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input}'


rule getChromSizes:
    input:
        rules.indexGenome.output.fai
    output:
        f'genome/{BUILD}-chrom.sizes'
    log:
        f'logs/getChromSizes/{BUILD}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


if config['ctcf'] is not None:

    rule getOrientationCTCF:
        input:
            genome = rules.bgzipGenome.output,
            index = rules.indexGenome.output,
            ctcf = lambda wc: config['ctcf'][int(wc.rep)]
        output:
            'genome/tracks/CTCF-{rep}.bed'
        log:
            'logs/getOrientationCTCF/{rep}.log'
        conda:
            f'{ENVS}/extract_ctcf_direction.yaml'
        shell:
            '{SCRIPTS}/extract_direction.py {input.genome} {input.ctcf} '
            '> {output} 2> {log}'


    rule processCTCFprediction:
        input:
            '/media/stephen/Data/genomes/hg19/tracks/GM12878-CTCF-hg19-prediction.tsv'
        output:
            #'genome/tracks/CTCF-{rep}.bed'
            'genome/tracks/CTCF-prediction.bed'
        params:
            threshold = 15
        log:
            'logs/processCTCFprediction.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/processCTCFBSDB.py --threshold {params.threshold} '
            '{input} > {output} 2> {log}'


    rule sortBed:
        input:
            #expand('genome/tracks/CTCF-{rep}.bed', rep = range(0, len(config['ctcf'])))
            rules.processCTCFprediction.output
        output:
            'genome/tracks/CTCF.sort.bed'
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
            'genome/tracks/CTCF.merged.bed'
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
            'genome/tracks/CTCF.scaled.bed'
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
            'genome/replicates/{rep}/tracks/CTCF-sampled-{rep}.bed'
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
            forward = 'genome/replicates/{rep}/tracks/CTCF-forward-{rep}.bed',
            reversed = 'genome/replicates/{rep}/tracks/CTCF-reverse-{rep}.bed'
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
            '--forward {output.forward} --reverse {output.reversed} '
            '{input} &> {log}'


rule scaleTracks:
    input:
        lambda wc: track_data[wc.track]['source']
    output:
        'genome/tracks/scaled/{track}'
    log:
        'logs/scaleTracks/{track}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/scaleBedScore.py --log {input} > {output} 2> {log}'


rule filterTracks:
    input:
        rules.scaleTracks.output
    output:
        'genome/replicates/{rep}/tracks/scaled/{track}'
    params:
        rep = REPS,
        seed = lambda wildcards: SEEDS[REPS.index(int(wildcards.rep))]
    log:
        'logs/filterTracks/{track}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/filterBedScore.py --seed {params.seed} '
        '{input} > {output} 2> {log}'


def getMasking(wc):
    """ Build masking command for maskFasta for track data """
    command = ''
    for track in track_data:
        character = track_data[track]['character']
        command += f'--bed genome/replicates/{wc.rep}/tracks/scaled/{track},{character} '
    if config['ctcf']:
        command += (f'--bed genome/replicates/{wc.rep}/tracks/CTCF-forward-{wc.rep}.bed,F '
                    f'--bed genome/replicates/{wc.rep}/tracks/CTCF-reverse-{wc.rep}.bed,R ')
    return command

def getMasking2(wc):
    """ Build masking command for maskFasta for track data """
    command = ''
    for bed, character in config['masking'].items():
        command += f'--bed {bed},{character} '
    if config['ctcf']:
        command += (f'--bed genome/replicates/{wc.rep}/tracks/CTCF-forward-{wc.rep}.bed,F '
                    f'--bed genome/replicates/{wc.rep}/tracks/CTCF-reverse-{wc.rep}.bed,R ')
    return command

rule maskFasta:
    input:
        rules.splitOrientation.output,
        expand('genome/replicates/{{rep}}/tracks/scaled/{track}',
            track = track_data.keys()),
        chromSizes = rules.getChromSizes.output
    output:
        pipe(f'genome/replicates/{{rep}}/sequence/{BUILD}-masked.fa')
    params:
        masking = getMasking
    group:
        'mask'
    log:
        f'logs/maskFasta/{BUILD}-{{rep}}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        '{SCRIPTS}/maskFastaV2.py {input.chromSizes} {params.masking} '
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
        f'logs/indexMasked/{BUILD}-{{rep}}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'


rule extractRegion:
    input:
        fasta = rules.bgzipMasked.output,
        index = rules.indexMasked.output
    output:
        pipe('genome/replicates/{rep}/genome/{region}-masked-{rep}.fa')
    group:
        'convert2Beads'
    params:
        region = f'{CHR}:{START}-{END}'
    log:
        'logs/extractRegion/{region}-{rep}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.fasta} {params.region} > {output} 2> {log}'


rule FastaToBeads:
    input:
        rules.extractRegion.output
    output:
        '{region}/{nbases}/reps/{rep}/sequence/{region}-beads-{rep}.txt'
    group:
        'convert2Beads'
    params:
        bases_per_bead = config['bases_per_bead']
    log:
        'logs/sequenceToBeads/{region}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/FastaToBeads.py --nbases {params.bases_per_bead} '
        '{input} > {output} 2> {log}'


rule BeadsToLammps:
    input:
        rules.FastaToBeads.output
    output:
        dat = '{region}/{nbases}/reps/{rep}/lammps/config/lammps_input.dat',
        coeffs = '{region}/{nbases}/reps/{rep}/lammps/config/coeffs.txt',
        groups = '{region}/{nbases}/reps/{rep}/lammps/config/groups.txt'
    params:
        nmonomers = config['monomers'],
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
        'logs/BeadsToLammps/{region}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_polymer.py --seed {params.seed} '
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
        '{region}/{nbases}/reps/{rep}/lammps/config/run-ctcf.lam'
    log:
        'logs/addCTCF/{region}-{nbases}-{rep}.log'
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
           '-var timestep {params.timestep} '
           '-var cosine_potential {params.cosine_potential} '
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
        warm_up = '{region}/{nbases}/reps/{rep}/lammps/warm_up.xyz.gz',
        simulation = '{region}/{nbases}/reps/{rep}/lammps/simulation.xyz.gz',
        radius_gyration = '{region}/{nbases}/reps/{rep}/lammps/radius_of_gyration.txt',
        restart = directory('{region}/{nbases}/reps/{rep}/lammps/restart/')
    params:
        timestep = config['timestep'],
        warm_up_time = config['warm_up'],
        sim_time = config['sim_time'],
        restart_time = int(config['timestep']) * 100, # Use 0 for no restart file
        cosine_potential = lambda wc: 10000 / config['bases_per_bead'],
        restart = lambda wc: f'{wc.region}/{wc.nbases}/reps/{wc.rep}/lammps/restart/Restart',
        seed = lambda wc: SEEDS[REPS.index(int(wc.rep))]
    group:
        'lammps'
    threads:
        config['threads'] if config['cluster'] else 1
    log:
        'logs/lammps/{region}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    shell:
        lmp_cmd


rule extractDNA:
    input:
        simulation = rules.lammps.output.simulation,
        groups = rules.BeadsToLammps.output.groups
    output:
        pipe('{region}/{nbases}/reps/{rep}/lammps/simulation-DNA.xyz')
    params:
        group = 'DNA'
    log:
        'logs/extractDNA/{region}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/filterXYZ.py --groupFile {input.groups} '
        '--group {params.group} <(zcat -f {input}) '
        '> {output} 2> {log}'



rule create_contact_matrix:
    input:
        rules.extractDNA.output
    output:
        '{region}/{nbases}/reps/{rep}/matrices/contacts.npz'
    params:
        distance =  3
    group:
        'lammps'
    log:
        'logs/create_contact_matrix/{region}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/create_contact_matrix.py '
        '--outdata {output} --distance {params.distance} '
        '{input} &> {log}'


rule mergeReplicates:
    input:
        expand('{{region}}/{{nbases}}/reps/{rep}/matrices/contacts.npz',
            rep=REPS)
    output:
        '{region}/{nbases}/merged/contacts.npz'
    params:
        method = config['method']
    log:
        'logs/mergeReplicates/{region}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/average_matrices.py --out {output} '
        '--method {params.method} {input} &> {log}'


rule matrix2homer:
    input:
        rules.mergeReplicates.output
    output:
        '{region}/{nbases}/merged/contacts.homer'
    params:
        chr = CHR,
        start = START,
        binsize = config['bases_per_bead']
    log:
        'logs/matrix2homer/{region}-{nbases}.log'
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
        '{region}/{nbases}/merged/contacts.h5'
    log:
        'logs/homer2H5/{region}-{nbases}.log'
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
        '{region}/{nbases}/merged/contacts-{binsize}.h5'
    params:
        nbins = MERGEBINS
    log:
        'logs/mergeBins/{region}-{nbases}-{binsize}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log}'


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
        '{region}/{nbases}/merged/configs-{binsize}.ini'
    conda:
        f'{ENVS}/python3.yaml'
    params:
        depth = int((END - START) / 2),
        hicConfig = getHiCconfig,
        colourMap = 'Purples'
    log:
        'logs/createConfig/{region}-{nbases}-{binsize}.log'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.matrix} '
        '--ctcf_orientation {input.ctcf_orientation} --log '
        '--colourmap {params.colourMap} '
        '--genes {input.genes} {params.hicConfig} '
        '--depth {params.depth} > {output} 2> {log}'


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        '{region}/{nbases}/merged/simulation-{binsize}.png'
    params:
        region = f'{CHR}:{START}-{END}',
        title = f'"{REGION} : {CHR}:{START}-{END}"',
        dpi = 600
    log:
        'logs/plotHiC/{region}-{nbases}-{binsize}.log'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    shell:
        'pyGenomeTracks --tracks {input} '
        '--region {params.region} '
        '--outFileName {output} '
        '--title {params.title} '
        '--dpi {params.dpi} &> {log}'


rule matrix2pre:
    input:
        rules.mergeReplicates.output
    output:
        '{region}/{nbases}/merged/contacts.pre.tsv'
    params:
        chr = CHR,
        start = START,
        binsize = config['bases_per_bead']
    log:
        'logs/matrix2pre/{region}-{nbases}.log'
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
        '{region}/{nbases}/merged/contacts.hic'
    log:
        'logs/juicerPre/{region}-{nbases}.log'
    params:
        chr = CHR,
        resolutions = '1000'
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

rule smooth_xyz:
    input:
        warm_up = '{region}/{nbases}/reps/1/lammps/warm_up.xyz.gz',
        simulation = '{region}/{nbases}/reps/1/lammps/simulation.xyz.gz'
    output:
        '{region}/{nbases}/reps/1/lammps/complete-smooth.xyz'
    params:
        window_size = config['window_size'],
        overlap = config['overlap']
    log:
        'logs/smooth_xyz/{region}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/smooth_xyz.py '
        '--size {params.window_size} '
        '--overlap {params.overlap} '
        '<(zcat {input.warm_up} {input.simulation}) '
        '> {output} 2> {log}'


checkpoint vmd:
    input:
        rules.smooth_xyz.output
    output:
        directory('{region}/{nbases}/vmd/sequence')
    log:
        'logs/vmd/{region}-{nbases}.log'
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
        '{region}/{nbases}/vmd/simulation.gif'
    params:
        delay = config['delay'],
        loop = config['loop']
    log:
        'logs/create_gif/{region}-{nbases}.log'
    conda:
        f'{ENVS}/imagemagick.yaml'
    shell:
        'convert -delay {params.delay} -loop {params.loop} '
        '{input.images} {output} &> {log}'
