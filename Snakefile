#!/usr/bin/env python3

container: "docker://continuumio/miniconda3:4.7.12"

import os
import random
import tempfile
from set_config import set_config, read_paths, nBeads

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
    'threads':        1,
    'ctcf':           None,
    'computeDirection': True,
    'masking':        {},
    'name' :          '',
    'genome' :        {'build':    'genome',
                       'sequence':  None,
                       'genes':     None,
                       'chr':       None,
                       'start':     None,
                       'end':       None,},
    'syntheticSequence' : None,
    'randomWalk':     False,
    'ignoreZeroPair': True,
    'min_rep':        1,
    'logScore':       True,
    'bases_per_bead': 1000,
    'monomers':       100,
    'method':         'mean',
    'n_molecules':    1000,
    'reps':           5,
    'random':         {'seed':             42,
                       'sequence':         True ,
                       'initialConform':   True ,
                       'monomerPositions': True,
                       'simulation':       True ,},
    'box':            {'xlo':        -50,
                       'xhi':         50,
                       'ylo':        -50,
                       'yhi':         50,
                       'zlo':        -50,
                       'zhi':         50,},
    'lammps':         {'restart':    0       ,
                       'timestep':   1000    ,
                       'warm_up':    20000   ,
                       'sim_time':   2000000 ,},
    'HiC':            {'matrix' :    None    ,
                       'binsize':    None    ,
                       'log' :       True    ,
                       'colourMap': 'Purples',
                       'dpi':        300 ,
                       'vMin':       None    ,
                       'vMax':       None    ,},
    'plotRG':         {'dpi':        300     ,
                       'confidence': 0.95    ,},
    'plotTU':         {'pvalue':     0.1     ,
                       'vMin':      -0.1     ,
                       'vMax':       0.1     ,},
    'delay':          10,
    'loop':           0,
    'tmpdir':         tempfile.gettempdir(),
    'coeffs':         '',
    'groupJobs':      False,
}
config = set_config(config, default_config)

workdir : config['workdir']

if config['syntheticSequence'] is None:
    # If no synthetic sequence ensure genome information is provided.
    invalid = False
    for key in ['sequence', 'chr', 'start', 'end']:
        if config['genome'][key] is None:
            sys.stderr.write(
                f'\033[31mNo synthetic sequence provided and '
                f'config["genome"]["{key}"] not set.\033[m\n')
            invalid = True
    if invalid:
        sys.exit('\033[31mInvalid configuration setting.\033[m\n')
else:
    config['genome']['sequence'] = []
    config['genome']['chr'] = 'synthetic'
    config['genome']['start'] = 1
    config['genome']['end'] = nBeads(config['syntheticSequence']) * config['bases_per_bead']

BUILD = config['genome']['build']
NAME = config['name']
CHR = config['genome']['chr']
START = config['genome']['start']
END = config['genome']['end']


# Define list of N reps from 1 to N
REPS = list(range(1, config['reps'] + 1))

# Set pipeline seed
random.seed(config['random']['seed'])
# Set seeds for generating bead sequence
if config['random']['sequence']:
    sequenceSeeds = [random.randint(1, (2**16) - 1) for rep in REPS]
else:
    sequenceSeeds = [random.randint(1, (2**16) - 1)] * config['reps']
# Set seeds for generating bead sequence
if config['random']['initialConform']:
    initialConformSeeds = [random.randint(1, (2**16) - 1) for rep in REPS]
else:
    initialConformSeeds = [random.randint(1, (2**16) - 1)] * config['reps']
# Set seeds for generating monomer positions
if config['random']['monomerPositions']:
    monomerSeeds = [random.randint(1, (2**16) - 1) for rep in REPS]
else:
    monomerSeeds = [random.randint(1, (2**16) - 1)] * config['reps']
# Set seeds for running lammps simulation
if config['random']['simulation']:
    simulationSeeds = [random.randint(1, (2**16) - 1) for rep in REPS]
else:
    simulationSeeds = [random.randint(1, (2**16) - 1)] * config['reps']


if config['HiC']['binsize'] is not None:
    BINSIZE = int(config['HiC']['binsize'])
    MERGEBINS, REMAINDER = divmod(BINSIZE, config['bases_per_bead'])
    if REMAINDER:
        sys.exit(
            f'\033[31Binsize {config["HiC"]["binsize"]} is not divisible by '
            f'bases per bead {config["bases_per_bead"]}.\033[m\n')
else:
    MERGEBINS = 1
    BINSIZE = config['bases_per_bead']


# Read track file basename and masking character into a dictionary
track_data = {}
for file, character in config['masking'].items():
    track_file = f'{os.path.basename(file)}'
    track_data[track_file] = {'source' : file, 'character' : character}

rule all:
    input:
        [expand('{name}/{nbases}/vmd/simulation.gif',
            nbases=config['bases_per_bead'], name=NAME),
         expand('{name}/{nbases}/merged/simulation-{binsize}.png',
            nbases=config['bases_per_bead'], name=NAME, binsize=BINSIZE),
         expand('{name}/{nbases}/merged/TU-{plot}.png',
            nbases=config['bases_per_bead'], name=NAME,
            plot=['Correlation','CircosPlot']),
         f'{NAME}/{config["bases_per_bead"]}/merged/radius_of_gyration.png']


rule unzipGenome:
    input:
        config['genome']['sequence']
    output:
        temp(f'genome/{BUILD}.fa')
    log:
        f'logs/unzipGenome/{BUILD}.log'
    shell:
        'zcat -f {input} > {output} 2> {log}'


rule indexGenome:
    input:
        rules.unzipGenome.output
    output:
        f'{rules.unzipGenome.output}.fai'
    log:
        f'logs/indexGenome/{BUILD}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input} &> {log}'


rule getChromSizes:
    input:
        rules.indexGenome.output
    output:
        f'genome/{BUILD}-chrom.sizes'
    log:
        f'logs/getChromSizes/{BUILD}.log'
    shell:
        'cut -f 1,2 {input} > {output} 2> {log}'


if config['ctcf'] is not None:

    if config['computeDirection']:

        rule modifyBedName:
            input:
                config['ctcf']
            output:
                'genome/tracks/CTCF-modifyName.bed'
            log:
                'logs/modifyBedName.log'
            conda:
                f'{ENVS}/python3.yaml'
            shell:
                '{SCRIPTS}/modifyBedName.py {input} > {output} 2> {log}'


        rule bed2Fasta:
            input:
                bed = rules.modifyBedName.output,
                genome = rules.unzipGenome.output,
                genomeIndex = rules.indexGenome.output
            output:
                'genome/tracks/CTCF-modifyName.fasta'
            log:
                'logs/bed2Fasta.log'
            conda:
                f'{ENVS}/bedtools.yaml'
            shell:
                'bedtools getfasta -name -fullHeader '
                '-bed {input.bed} -fi {input.genome} > {output} 2> {log}'


        rule runCTCFpredict:
            input:
                rules.bed2Fasta.output
            output:
                'genome/tracks/CTCF-prediction.tsv'
            log:
                'logs/runCTCFpredict.log'
            conda:
                f'{ENVS}/python3.yaml'
            shell:
                '{SCRIPTS}/runCTCFpredict.py {input} > {output} 2> {log}'


        rule processCTCFprediction:
            input:
                rules.runCTCFpredict.output
            output:
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

    def CTCFinput(wildcards):
        if config['computeDirection']:
            return rules.processCTCFprediction.output
        else:
            return config['ctcf']

    rule sortBed:
        input:
            CTCFinput
        output:
            'genome/tracks/CTCF.sort.bed'
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
        params:
            logScore = '--log' if config['logScore'] else ''
        group:
            'bedtools'
        log:
            'logs/scaleBed.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/scaleBedScore.py {params.logScore} {input} '
            '> {output} 2> {log}'


    rule filterBedScore:
        input:
            rules.scaleBed.output
        output:
            'genome/replicates/{rep}/tracks/CTCF-sampled-{rep}.bed'
        params:
            rep = REPS,
            seed = lambda wc: sequenceSeeds[int(wc.rep) - 1]
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
    params:
        logScore = '--log' if config['logScore'] else ''
    log:
        'logs/scaleTracks/{track}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/scaleBedScore.py {params.logScore} {input} '
        '{params.logScore} > {output} 2> {log}'


rule filterTracks:
    input:
        rules.scaleTracks.output
    output:
        'genome/replicates/{rep}/tracks/scaled/{track}'
    params:
        rep = REPS,
        seed = lambda wc: sequenceSeeds[int(wc.rep) - 1]
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


def getSplitOrient(wc):
    """ Return rule output only if used. """
    if config['ctcf']:
        return rules.splitOrientation.output
    else:
        return []

rule maskFasta:
    input:
        getSplitOrient,
        expand('genome/replicates/{{rep}}/tracks/scaled/{track}',
            track = track_data.keys()),
        chromSizes = rules.getChromSizes.output
    output:
        pipe(f'genome/replicates/{{rep}}/sequence/{BUILD}-masked.fa')
    params:
        masking = getMasking,
        seed = lambda wc: sequenceSeeds[int(wc.rep) - 1]
    group:
        'mask'
    log:
        f'logs/maskFasta/{BUILD}-{{rep}}.log'
    conda:
        f'{ENVS}/bedtools.yaml'
    shell:
        '{SCRIPTS}/maskFastaV2.py {input.chromSizes} {params.masking} '
        '--seed {params.seed} > {output} 2> {log}'


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
        pipe('genome/replicates/{rep}/genome/{name}-masked-{rep}.fa')
    group:
        'convert2Lammps'
    params:
        region = f'{CHR}:{START}-{END}'
    log:
        'logs/extractRegion/{name}-{rep}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        'samtools faidx {input.fasta} {params.region} > {output} 2> {log}'


rule FastaToBeads:
    input:
        rules.extractRegion.output
    output:
        '{name}/{nbases}/reps/{rep}/sequence/{name}-beads-{rep}.txt'
    group:
        'convert2Lammps'
    params:
        bases_per_bead = config['bases_per_bead'],
        seed = lambda wc: sequenceSeeds[int(wc.rep) - 1]
    log:
        'logs/sequenceToBeads/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/FastaToBeads.py --nbases {params.bases_per_bead} '
        '--seed {params.seed} {input} > {output} 2> {log}'


def beadsInput(wc):
    if config['syntheticSequence'] is not None:
        return config['syntheticSequence']
    else:
        return rules.FastaToBeads.output

rule BeadsToLammps:
    input:
        beadsInput
    output:
        dat = '{name}/{nbases}/reps/{rep}/lammps/config/lammps_input.dat',
        coeffs = '{name}/{nbases}/reps/{rep}/lammps/config/coeffs.txt',
        groups = '{name}/{nbases}/reps/{rep}/lammps/config/groups.txt'
    params:
        nMonomers = config['monomers'],
        coeffs = config['coeffs'],
        basesPerBead = config['bases_per_bead'],
        polymerSeed = lambda wc: initialConformSeeds[int(wc.rep) - 1],
        monomerSeed = lambda wc: monomerSeeds[int(wc.rep) - 1],
        xlo = config['box']['xlo'],
        xhi = config['box']['xhi'],
        ylo = config['box']['ylo'],
        yhi = config['box']['yhi'],
        zlo = config['box']['zlo'],
        zhi = config['box']['zhi'],
        randomWalk = '--randomWalk' if config['randomWalk'] else ''
    group:
        'lammps'
    log:
        'logs/BeadsToLammps/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_polymer.py '
        '--polymerSeed {params.polymerSeed} '
        '--monomerSeed {params.monomerSeed} '
        '--xlo {params.xlo} --xhi {params.xhi} '
        '--ylo {params.ylo} --yhi {params.yhi} '
        '--zlo {params.zlo} --zhi {params.zhi} '
        '--ctcf --coeffOut {output.coeffs} '
        '--groupOut {output.groups} '
        '--nMonomer {params.nMonomers} '
        '--pairCoeffs {params.coeffs} '
        '--basesPerBead {params.basesPerBead} '
        '{params.randomWalk} {input} > {output.dat} 2> {log}'


rule checkCTCF:
    input:
        script = f'{SCRIPTS}/run-ctcf.lam',
        dat = rules.BeadsToLammps.output.dat
    output:
        temp('{name}/{nbases}/reps/{rep}/lammps/config/run-ctcf.tmp.lam')
    group:
        'lammps'
    log:
        'logs/checkCTCF/{name}-{nbases}-{rep}.log'
    shell:
        "(awk '/Bonds/,/Angles/{{print $2}}' {input.dat} "
        "| grep -q 2 && cat {input.script} > {output} "
        "|| sed -e '/# CTCF/d' {input.script} "
        "> {output}) 2> {log}"


rule addCTCF:
    input:
        script = rules.checkCTCF.output,
        coeffs = rules.BeadsToLammps.output.coeffs,
        groups = rules.BeadsToLammps.output.groups,
    output:
        '{name}/{nbases}/reps/{rep}/lammps/config/run-ctcf.lam'
    group:
        'lammps'
    log:
        'logs/addCTCF/{name}-{nbases}-{rep}.log'
    shell:
        "sed '/^#PAIR_COEFF/ r {input.coeffs}' {input.script} "
        "| sed '/^#GROUPS/ r {input.groups}' > {output} 2> {log}"


lmp_cmd = ('-var infile {input.data} '
           '-var radius_gyration {output.radius_gyration} '
           '-var warm_up {output.warm_up} '
           '-var warm_up_time {params.warm_up_time} '
           '-var restart {params.restart} '
           '-var sim {output.simulation} '
           '-var sim_time {params.sim_time} '
           '-var timestep {params.timestep} '
           '-var cosine_potential {params.cosine_potential} '
           '-var seed {params.seed} '
           '-in {input.script} '
           '-log /dev/null &> {log}')
if config['cluster']:
    #lmp_cmd = 'mpirun -n {threads} -x OMP_NUM_THREADS=1 lmp_mpi ' + lmp_cmd
    lmp_cmd = 'mpirun -np {threads} lmp_mpi ' + lmp_cmd
else:
    lmp_cmd = 'lmp_serial ' + lmp_cmd


def restartCommand(wc):
    # Define reset command, turn restart OFF if set to 0
    step = int(config['lammps']['restart'])
    if step == 0:
        return f'{step}'
    else:
        prefix = f'{wc.region}/{wc.nbases}/reps/{wc.rep}/lammps/restart/Restart'
        interval = step * int(config['lammps']['timestep'])
        return f'"{interval} {prefix}"'

rule lammps:
    input:
        data = rules.BeadsToLammps.output.dat,
        script = rules.addCTCF.output
    output:
        warm_up = '{name}/{nbases}/reps/{rep}/lammps/warm_up.custom.gz',
        simulation = '{name}/{nbases}/reps/{rep}/lammps/simulation.custom.gz',
        radius_gyration = '{name}/{nbases}/reps/{rep}/lammps/radius_of_gyration.txt',
        restart = directory('{name}/{nbases}/reps/{rep}/lammps/restart/')
    params:
        timestep = config['lammps']['timestep'],
        warm_up_time = config['lammps']['warm_up'],
        sim_time = config['lammps']['sim_time'],
        cosine_potential = lambda wc: 10000 / config['bases_per_bead'],
        restart = restartCommand,
        seed = lambda wc: simulationSeeds[int(wc.rep) - 1],
    group:
        'lammps'
    threads:
        config['threads'] if config['cluster'] else 1
    log:
        'logs/lammps/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    shell:
        lmp_cmd


rule plotRG:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/lammps/radius_of_gyration.txt',
            rep=REPS)
    output:
        '{name}/{nbases}/merged/radius_of_gyration.png'
    params:
        confidence = config['plotRG']['confidence'],
        dpi = config['plotRG']['dpi']
    group:
        'lammps' if config['groupJobs'] else 'plotRG'
    log:
        'logs/plotRG/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotRG.py --output {output} --dpi {params.dpi} '
        '--confidence {params.confidence} {input} 2> {log}'


rule custom2XYZ:
    input:
        '{name}/{nbases}/reps/{rep}/lammps/{mode}.custom.gz'
    output:
        temp('{name}/{nbases}/reps/{rep}/lammps/{mode}.xyz.gz')
    group:
        'processAllLammps' if config['groupJobs'] else 'custom2XYZ'
    log:
        'logs/custom2XYZ/{name}-{nbases}-{rep}-{mode}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '(zcat {input} | {SCRIPTS}/custom2XYZ.py | gzip > {output}) 2> {log}'


rule getAtomGroups:
    input:
        rules.BeadsToLammps.output.dat
    output:
        '{name}/{nbases}/reps/{rep}/lammps/config/atomGroups.json'
    group:
        'processAllLammps' if config['groupJobs'] else 'getAtomGroups'
    log:
        'logs/getAtomGroups/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/getAtomGroups.py {input} > {output} 2> {log}'


rule filterXYZ:
    input:
        simulation = '{name}/{nbases}/reps/{rep}/lammps/simulation.xyz.gz',
        groups = rules.getAtomGroups.output
    output:
        '{name}/{nbases}/reps/{rep}/lammps/simulation-{bead}.xyz.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'filterXYZ'
    log:
        'logs/filterXYZ/{name}-{nbases}-{rep}-{bead}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/filterXYZ.py {wildcards.bead} '
        '{input.groups} <(zcat -f {input.simulation}) '
        '| gzip > {output} 2> {log}'


rule computeTUCorrelation:
    input:
        monomers = '{name}/{nbases}/reps/{rep}/lammps/simulation-TF.xyz.gz',
        dna = '{name}/{nbases}/reps/{rep}/lammps/simulation-TU.xyz.gz',
        groups = rules.getAtomGroups.output
    output:
        '{name}/{nbases}/reps/{rep}/TUcorrelation.csv'
    params:
        distance = 1.8,
        ignoreZeroPair = '--ignoreZeroPair' if config['ignoreZeroPair'] else ''
    group:
        'processAllLammps' if config['groupJobs'] else 'computeTUCorrelation'
    log:
        'logs/computeTUCorrelation/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/computeTUCorrelation.py '
        '<(zcat -f {input.dna}) <(zcat -f {input.monomers}) {input.groups} '
        '--distance {params.distance} {params.ignoreZeroPair} '
        '--out {output} &> {log}'


rule plotTUCorrelation:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/TUcorrelation.csv', rep=REPS),
    output:
        heatmap = '{name}/{nbases}/merged/TU-Correlation.png',
        circos = '{name}/{nbases}/merged/TU-CircosPlot.png'
    params:
        pvalue = config['plotTU']['pvalue'],
        vMin = config['plotTU']['vMin'],
        vMax = config['plotTU']['vMax']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotTUCorrelation'
    log:
        'logs/plotTUCorrelation/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotTUCorrelation.py {input} --heatmap {output.heatmap} '
        '--circos {output.circos} --pvalue {params.pvalue} '
        '--vmin {params.vMin} --vmax {params.vMax} &> {log}'


rule createContactMatrix:
    input:
        '{name}/{nbases}/reps/{rep}/lammps/simulation-DNA.xyz.gz'
    output:
        '{name}/{nbases}/reps/{rep}/matrices/contacts.npz'
    params:
        distance =  3
    group:
        'processAllLammps' if config['groupJobs'] else 'createContactMatrix'
    log:
        'logs/createContactMatrix/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/create_contact_matrix.py '
        '--outdata {output} --distance {params.distance} '
        '<(zcat -f {input}) &> {log}'


rule mergeReplicates:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/matrices/contacts.npz',
            rep=REPS)
    output:
        '{name}/{nbases}/merged/contacts.npz'
    params:
        method = config['method']
    group:
        'processAllLammps' if config['groupJobs'] else 'createContactMatrix'
    log:
        'logs/mergeReplicates/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/average_matrices.py --out {output} '
        '--method {params.method} {input} &> {log}'


rule matrix2homer:
    input:
        rules.mergeReplicates.output
    output:
        '{name}/{nbases}/merged/contacts.homer'
    params:
        chr = CHR,
        start = START,
        binsize = config['bases_per_bead']
    log:
        'logs/matrix2homer/{name}-{nbases}.log'
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
        '{name}/{nbases}/merged/contacts.h5'
    log:
        'logs/homer2H5/{name}-{nbases}.log'
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
        '{name}/{nbases}/merged/contacts-{binsize}.h5'
    params:
        nbins = MERGEBINS
    log:
        'logs/mergeBins/{name}-{nbases}-{binsize}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicMergeMatrixBins --matrix {input} --numBins {params.nbins} '
        '--outFileName {output} &> {log}'


def getHiCconfig(wc):
    """ Build config command for configuration """
    command = ''
    if not config['syntheticSequence']:
        if config['HiC']['matrix']:
            command += f'--flip --matrix2 {config["HiC"]["matrix"]} '
        if config['genome']['genes']:
            command += f' --genes {config["genome"]["genes"]}'
    if config['HiC']['vMin'] is not None:
        command += f' --vMin {config["plot"]["vMin"]} '
    if config['HiC']['vMax'] is not None:
        command += f' --vMax {config["plot"]["vMax"]} '
    if config['HiC']['log']:
        command += ' --log '
    if config['syntheticSequence'] is not None and config['ctcf'] is not None:
        command += f'--ctcfOrient {rules.scaleBed.output}'

    return command


def getCTCFOrient(wc):
    """ Return CTCF orientation if not synthetic """
    if config['syntheticSequence'] is None:
        return rules.scaleBed.output
    else:
        return []

rule createConfig:
    input:
        matrix = rules.mergeBins.output,
        ctcfOrient = getCTCFOrient
    output:
        '{name}/{nbases}/merged/configs-{binsize}.ini'
    params:
        depth = END - START + 1,
        hicConfig = getHiCconfig,
        colourMap = config['HiC']['colourMap']
    group:
        'plotHiC'
    log:
        'logs/createConfig/{name}-{nbases}-{binsize}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.matrix} '
        '--colourMap {params.colourMap} --depth {params.depth} '
        '{params.hicConfig} > {output} 2> {log}'


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        '{name}/{nbases}/merged/simulation-{binsize}.png'
    params:
        region = f'{CHR}:{START}-{END}',
        title = f'"{NAME} : {CHR}:{START}-{END}"',
        dpi = config['HiC']['dpi']
    group:
        'plotHiC'
    log:
        'logs/plotHiC/{name}-{nbases}-{binsize}.log'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    shell:
        'pyGenomeTracks --tracks {input} '
        '--region {params.region} '
        '--outFileName {output} '
        '--title {params.title} '
        '--dpi {params.dpi} &> {log}'


if config['syntheticSequence'] is None:

    rule matrix2pre:
        input:
            rules.mergeReplicates.output
        output:
            '{name}/{nbases}/merged/contacts.pre.tsv'
        params:
            chr = CHR,
            start = START,
            binsize = config['bases_per_bead']
        log:
            'logs/matrix2pre/{name}-{nbases}.log'
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
            '{name}/{nbases}/merged/contacts.hic'
        log:
            'logs/juicerPre/{name}-{nbases}.log'
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


checkpoint vmd:
    input:
        '{name}/{nbases}/reps/1/lammps/simulation.xyz.gz'
    output:
        directory('{name}/{nbases}/vmd/sequence')
    log:
        'logs/vmd/{name}-{nbases}.log'
    conda:
        f'{ENVS}/vmd.yaml'
    shell:
        'vmd -eofexit -nt -dispdev text -args <(zcat {input}) {output} '
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
        '{name}/{nbases}/vmd/simulation.gif'
    params:
        delay = config['delay'],
        loop = config['loop']
    log:
        'logs/create_gif/{name}-{nbases}.log'
    conda:
        f'{ENVS}/imagemagick.yaml'
    shell:
        'convert -delay {params.delay} -loop {params.loop} '
        '{input.images} {output} &> {log}'
