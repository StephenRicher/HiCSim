#!/usr/bin/env python3

container: "docker://continuumio/miniconda3:4.7.12"

import os
import random
import tempfile
from set_config import set_config, read_paths, nBeads, adjustCoordinates

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
    'workdir':        workflow.basedir,
    'ctcf':           {'data':             None,
                       'computeDirection': True},
    'atac':           {'bedgraph':    None  ,
                       'scale':      'none' ,
                       'percentile':  97.5  ,},
    'masking':        {},
    'genome' :        {'build':    'genome',
                       'sequence':  None,
                       'name':      None,
                       'genes':     None,
                       'chr':       None,
                       'start':     None,
                       'end':       None,},
    'syntheticSequence' : {}              ,
    'min_rep':        1,
    'scaleBed':      'sqrt',
    'bases_per_bead': 1000,
    'monomers':       100,
    'method':         'mean',
    'coeffs':         '',
    'reps':           5,
    'random':         {'seed':             42,
                       'walk':             False,
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
                       'dpi':        300     ,
                       'vMin':       None    ,
                       'vMax':       None    ,
                       'vMin2':      None    ,
                       'vMax2':      None    ,},
    'dtwDNA':         {'sampletime': 0       ,
                       'maxtime':    100000 ,},
    'dtwTU':          {'sampletime': 0       ,
                       'maxtime':    0       ,},
    'plotRG':         {'dpi':        300     ,
                       'confidence': 0.95    ,},
    'plotTU':         {'pvalue':     0.1     ,
                       'vMin':      -0.3     ,
                       'vMax':       0.3     ,
                       'minRep':     1       ,
                       'fontSize':   14      ,},
    'GIF':            {'create':     True    ,
                       'delay':      10      ,
                       'loop':       0       ,},
    'groupJobs':      False,
    'cluster':        False,
}
config = set_config(config, default_config)

workdir : config['workdir']

if config['scaleBed'] not in ['sqrt', 'none', 'log']:
    sys.exit('"scaleBed" in configuration file must be one of either '
             '"none", "sqrt" or "log".')

details = {}
if not config['syntheticSequence']:
    # If no synthetic sequence ensure genome information is provided.
    invalid = False
    for key in ['sequence', 'name', 'chr', 'start', 'end']:
        if config['genome'][key] is None:
            sys.stderr.write(
                f'\033[31mNo synthetic sequence provided and '
                f'config["genome"]["{key}"] not set.\033[m\n')
            invalid = True
    if invalid:
        sys.exit('\033[31mInvalid configuration setting.\033[m\n')
    for name in [config['genome']['name']]:
        start, end = adjustCoordinates(
            config['genome']['start'],
            config['genome']['end'],
            config['bases_per_bead'])
        scaledRegion = f"{config['genome']['chr']}-{start}-{end}"
        print(f'Adjusting {name} positions to {scaledRegion}')
        details[name] = {'chr':   config['genome']['chr'],
                         'start': start, 'end': end}
else:
    config['genome']['sequence'] = []
    for name in config['syntheticSequence'].keys():
        end = (nBeads(config['syntheticSequence'][name])
               * config['bases_per_bead'])
        details[name] = {'chr':   name,
                         'start': 1,
                         'end':   end}

BUILD = config['genome']['build']

# Define list of N reps from 1 to N
REPS = list(range(1, config['reps'] + 1))

# Set seeds for different parts of workflow.
seeds = {}
for type in ['sequence', 'initialConform', 'monomerPositions', 'simulation']:
    random.seed(config['random']['seed'])
    if config['random'][type]:
        seeds[type] = [random.randint(1, (2**16) - 1) for rep in REPS]
    else:
        seeds[type] = [random.randint(1, (2**16) - 1)] * config['reps']

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
if config['masking']:
    for file, character in config['masking'].items():
        track_file = f'{os.path.basename(file)}'
        track_data[track_file] = {'source' : file, 'character' : character}

wildcard_constraints:
    all = r'[^\/]+',
    name = rf'{"|".join(details.keys())}',
    binsize = rf'{BINSIZE}',
    nbases = rf'{config["bases_per_bead"]}',
    rep = rf'{"|".join([str(rep) for rep in REPS])}'


rule all:
    input:
        [expand('vmd/{name}-{nbases}-1-simulation.gif',
            nbases=config['bases_per_bead'],
            name=details.keys()) if config['GIF']['create'] else [],
         expand('plots/contactMatrix/{name}-{nbases}-{binsize}-contactMatrix.png',
            nbases=config['bases_per_bead'], name=details.keys(), binsize=BINSIZE),
         expand('plots/{plot}/{name}-{nbases}-{plot}.png',
            nbases=config['bases_per_bead'], name=details.keys(),
            plot=['pairCluster', 'meanVariance','TUcorrelation', 'TUcircos',
                  'TUactivation', 'radiusGyration', 'TUreplicateCount',])]


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


if config['ctcf']['data'] is not None:

    if config['ctcf']['computeDirection']:

        rule modifyBedName:
            input:
                config['ctcf']['data']
            output:
                'tracks/CTCF-modifyName.bed'
            group:
                'processCTCF'
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
                'tracks/CTCF-modifyName.fasta'
            group:
                'processCTCF'
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
                'tracks/CTCF-prediction.tsv'
            group:
                'processCTCF'
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
                'tracks/CTCF-prediction.bed'
            group:
                'processCTCF'
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
        if config['ctcf']['computeDirection']:
            return rules.processCTCFprediction.output
        else:
            return config['ctcf']['data']


    rule sortCTCF:
        input:
            CTCFinput
        output:
            'tracks/CTCF-sort.bed'
        group:
            'processCTCF'
        log:
            'logs/sortBed.log'
        conda:
            f'{ENVS}/bedtools.yaml'
        shell:
            'bedtools sort -i {input} > {output} 2> {log}'


    rule mergeCTCF:
        input:
            rules.sortCTCF.output
        output:
            'tracks/CTCF-merged.bed'
        group:
            'processCTCF'
        log:
            'logs/mergeBed.log'
        conda:
            f'{ENVS}/bedtools.yaml'
        shell:
            'bedtools merge -s -c 4,5,6 -o count,median,distinct '
            '-i {input} > {output} 2> {log}'


    rule scaleCTCF:
        input:
            rules.mergeCTCF.output
        output:
            'tracks/CTCF-scaled.bed'
        params:
            transform =  config['scaleBed']
        group:
            'processCTCF'
        log:
            'logs/scaleBed.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/scaleBedScore.py --transform {params.transform} '
            '{input} > {output} 2> {log}'


    rule filterCTCF:
        input:
            rules.scaleCTCF.output
        output:
            'tracks/{rep}/CTCF-sampled.bed'
        params:
            rep = REPS,
            seed = lambda wc: seeds['sequence'][int(wc.rep) - 1]
        group:
            'lammps'
        log:
            'logs/sampleBed/{rep}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '{SCRIPTS}/filterBedScore.py --seed {params.seed} '
            '{input} > {output} 2> {log}'


    rule splitOrientation:
        input:
            rules.filterCTCF.output
        output:
            forward = 'tracks/{rep}/CTCF-forward.bed',
            reversed = 'tracks/{rep}/CTCF-reverse.bed'
        params:
            min_rep = config['min_rep']
        group:
            'lammps'
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
        'tracks/scaled/{track}'
    params:
        transform = config['scaleBed']
    group:
        'lammps'
    log:
        'logs/scaleTracks/{track}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/scaleBedScore.py --transform {params.transform} '
        '{input} > {output} 2> {log}'

rule processATAC:
    input:
        config['atac']['bedgraph']
    output:
        'tracks/ATAC/ATAC-beadModifier-{nbases}.json'
    params:
        transform = config['atac']['scale'],
        percentile = config['atac']['percentile'],
        precision = 2
    group:
        'lammps'
    log:
        'logs/processATAC-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/processATAC.py --transform {params.transform} '
        '--nbases {wildcards.nbases} --percentile {params.percentile} '
        '--precision {params.precision} {input} > {output} 2> {log}'


def getRegion(wc):
    chrom = details[wc.name]['chr']
    start = details[wc.name]['start']
    end = details[wc.name]['end']
    return f'{chrom}:{start}-{end}'


rule subsetATAC:
    input:
        rules.processATAC.output
    output:
        '{name}/{nbases}/ATAC-beadModifier.json'
    params:
        region = getRegion
    group:
        'lammps'
    log:
        'logs/subsetATAC/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/subsetATAC.py --nbases {wildcards.nbases} '
        '--region {params.region} {input} > {output} 2> {log}'


rule filterTracks:
    input:
        rules.scaleTracks.output
    output:
        'tracks/scaled/{rep}/{track}'
    params:
        rep = REPS,
        seed = lambda wc: seeds['sequence'][int(wc.rep) - 1]
    group:
        'lammps'
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
        command += f'--bed tracks/scaled/{wc.rep}/{track},{character} '
    if config['ctcf']['data']:
        command += (f'--bed tracks/{wc.rep}/CTCF-forward.bed,F '
                    f'--bed tracks/{wc.rep}/CTCF-reverse.bed,R ')
    return command


def getSplitOrient(wc):
    """ Return rule output only if used. """
    if config['ctcf']['data']:
        return rules.splitOrientation.output
    else:
        return []


rule maskFasta:
    input:
        getSplitOrient,
        expand('tracks/scaled/{{rep}}/{track}', track=track_data.keys()),
        chromSizes = rules.getChromSizes.output
    output:
        pipe('{name}/reps/{name}-{rep}-masked.fa')
    params:
        masking = getMasking,
        region = getRegion,
        seed = lambda wc: seeds['sequence'][int(wc.rep) - 1]
    group:
        'lammps'
    log:
        'logs/maskFasta/maskFasta-{name}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/maskFastaV2.py {input.chromSizes} {params.masking} '
        '--region {params.region} --seed {params.seed} > {output} 2> {log}'


rule FastaToBeads:
    input:
        rules.maskFasta.output
    output:
        '{name}/{nbases}/reps/{rep}/{name}-beads-{rep}.txt'
    group:
        'lammps'
    params:
        bases_per_bead = config['bases_per_bead'],
        seed = lambda wc: seeds['sequence'][int(wc.rep) - 1]
    log:
        'logs/sequenceToBeads/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/FastaToBeads.py --nbases {params.bases_per_bead} '
        '--seed {params.seed} {input} > {output} 2> {log}'


def beadsInput(wc):
    if config['syntheticSequence']:
        return config['syntheticSequence'][wc.name]
    else:
        return rules.FastaToBeads.output


def setBeadsToLammpsCmd():
    cmd = '{SCRIPTS}/generate_polymer.py '
    if config['atac']['bedgraph']:
        cmd += '--atac {input.atac} '
    cmd += (
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
        '{params.randomWalk} {input.beads} > {output.dat} 2> {log}')
    return cmd


rule BeadsToLammps:
    input:
        beads = beadsInput,
        atac = rules.subsetATAC.output if config['atac']['bedgraph'] else []
    output:
        dat = '{name}/{nbases}/reps/{rep}/lammps/config/lammps_input.dat',
        coeffs = '{name}/{nbases}/reps/{rep}/lammps/config/coeffs.txt',
        groups = '{name}/{nbases}/reps/{rep}/lammps/config/groups.txt'
    params:
        nMonomers = config['monomers'],
        coeffs = config['coeffs'],
        basesPerBead = config['bases_per_bead'],
        polymerSeed = lambda wc: seeds['initialConform'][int(wc.rep) - 1],
        monomerSeed = lambda wc: seeds['monomerPositions'][int(wc.rep) - 1],
        xlo = config['box']['xlo'],
        xhi = config['box']['xhi'],
        ylo = config['box']['ylo'],
        yhi = config['box']['yhi'],
        zlo = config['box']['zlo'],
        zhi = config['box']['zhi'],
        randomWalk = '--randomWalk' if config['random']['walk'] else ''
    group:
        'lammps'
    log:
        'logs/BeadsToLammps/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        setBeadsToLammpsCmd()


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
if config['cluster'] and not config['groupJobs']:
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
        seed = lambda wc: seeds['simulation'][int(wc.rep) - 1],
    group:
        'lammps'
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
        'plots/radiusGyration/{name}-{nbases}-radiusGyration.png'
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


rule writeTUdistribution:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/lammps/config/atomGroups.json', rep=REPS)
    output:
        '{name}/{nbases}/merged/info/{name}-TU-distribution.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'writeTUdistribution'
    log:
        'logs/writeTUdistribution/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '({SCRIPTS}/writeTUdistribution.py {input} | gzip > {output}) &> {log}'


rule reformatLammps:
    input:
        rules.lammps.output.simulation
    output:
        '{name}/{nbases}/reps/{rep}/lammps/simulation.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'reformatLammps'
    log:
        'logs/reformatLammps/{name}-{nbases}-{rep}.log'
    shell:
        '({SCRIPTS}/reformatLammps.awk <(zcat {input}) | gzip > {output}) &> {log}'


rule processTUinfo:
    input:
        sim = rules.reformatLammps.output,
        groups = rules.getAtomGroups.output
    output:
        '{name}/{nbases}/reps/{rep}/TU-info.csv.gz',
    params:
        distance = 1.8
    group:
        'processAllLammps' if config['groupJobs'] else 'processTUinfo'
    log:
        'logs/processTUinfo/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/processTUinfo.py --out {output} '
        '--distance {params.distance} {input.groups} {input.sim} &> {log}'


rule DBSCAN:
    input:
        rules.processTUinfo.output
    output:
        clusterPairs = '{name}/{nbases}/reps/{rep}/TU-clusterPair.csv.gz',
        clusterPlot = 'plots/DBSCAN/{name}/{name}-{nbases}-{rep}-cluster.png',
    params:
        eps = 6,
        minSamples = 2
    group:
        'DBSCAN_all' if config['groupJobs'] else 'DBSCAN'
    log:
        'logs/DBSCAN/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/dbscanTU.py --outplot {output.clusterPlot} '
        '--eps {params.eps} --minSamples {params.minSamples} '
        '--out {output.clusterPairs} {input} &> {log}'


rule plotDBSCAN:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/TU-clusterPair.csv.gz', rep=REPS),
    output:
         'plots/pairCluster/{name}-{nbases}-pairCluster.png'
    group:
        'DBSCAN_all' if config['groupJobs'] else 'plotDBSCAN'
    log:
        'logs/plotDBSCAN/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotDBSCAN.py {output} {input} &> {log}'


rule computeTUstats:
    input:
        '{name}/{nbases}/reps/{rep}/TU-info.csv.gz'
    output:
        '{name}/{nbases}/reps/{rep}/TU-stats.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'computeTUstats'
    log:
        'logs/computeTUstats/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/computeTUstats.py --out {output} {input} &> {log}'


rule mergeByRep:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/TU-stats.csv.gz', rep=REPS),
    output:
        '{name}/{nbases}/merged/{name}-TU-stats.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'mergeByRep'
    log:
        'logs/mergeByRep/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mergeByRep.py --out {output} {input} &> {log}'


rule plotMeanVariance:
    input:
        rules.mergeByRep.output
    output:
        'plots/meanVariance/{name}-{nbases}-meanVariance.png'
    params:
        fontSize = config['plotTU']['fontSize']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotMeanVariance'
    log:
        'logs/plotMeanVariance/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotMeanVariance.py {input} {output} '
        '--fontSize {params.fontSize}  &> {log}'


rule computeTUcorrelation:
    input:
        rules.mergeByRep.output
    output:
        '{name}/{nbases}/merged/{name}-TU-correlation.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'computeTUcorrelation'
    log:
        'logs/computeTUcorrelation/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/computeTUcorrelation.py --out {output} {input} &> {log}'


rule plotTUcorrelation:
    input:
        correlations = rules.computeTUcorrelation.output,
        beadDistribution = rules.writeTUdistribution.output
    output:
        meanHeatmap = 'plots/TUcorrelation/{name}-{nbases}-TUcorrelation.png',
        sumHeatmap = 'plots/TUreplicateCount/{name}-{nbases}-TUreplicateCount.png',
        circos = 'plots/TUcircos/{name}-{nbases}-TUcircos.png'
    params:
        pvalue = config['plotTU']['pvalue'],
        vMin = config['plotTU']['vMin'],
        vMax = config['plotTU']['vMax'],
        minRep = config['plotTU']['minRep'],
        fontSize = config['plotTU']['fontSize']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotTUcorrelation'
    log:
        'logs/plotTUcorrelation/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotTUcorrelation.py {input.beadDistribution} '
        '{input.correlations} --circos {output.circos} '
        '--sumHeatmap {output.sumHeatmap} --meanHeatmap {output.meanHeatmap} '
        '--fontSize {params.fontSize} --pvalue {params.pvalue} '
        '--minRep {params.minRep} --vmin {params.vMin} '
        '--vmax {params.vMax} &> {log}'


rule plotTUactivation:
    input:
        expand(
            '{{name}}/{{nbases}}/reps/{rep}/TU-info.csv.gz', rep=REPS),
    output:
        'plots/TUactivation/{name}-{nbases}-TUactivation.png'
    group:
        'processAllLammps' if config['groupJobs'] else 'plotTUactivation'
    log:
        'logs/plotTUactivation/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotTUactivation.py {input} --out {output} &> {log}'


rule createContactMatrix:
    input:
        xyz = rules.reformatLammps.output,
        groups = rules.getAtomGroups.output
    output:
        '{name}/{nbases}/reps/{rep}/matrices/contacts.npz'
    params:
        distance = 3,
        periodic = '--periodic',
        x = abs(config['box']['xhi'] - config['box']['xlo']),
        y = abs(config['box']['yhi'] - config['box']['ylo']),
        z = abs(config['box']['zhi'] - config['box']['zlo'])
    group:
        'processAllLammps' if config['groupJobs'] else 'createContactMatrix'
    log:
        'logs/createContactMatrix/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/create_contact_matrix.py {params.periodic} '
        '--outdata {output} --distance {params.distance} '
        '--dimensions {params.x} {params.y} {params.z} {input.groups} '
        '<(zcat -f {input.xyz}) &> {log}'


rule mergeReplicates:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/matrices/contacts.npz',
            rep=REPS)
    output:
        '{name}/{nbases}/merged/matrices/{name}.npz'
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
        '{name}/{nbases}/merged/matrices/{name}.npz'
    output:
        '{name}/{nbases}/merged/matrices/{name}.homer'
    params:
        chr = lambda wc: details[wc.name]['chr'],
        start = lambda wc: details[wc.name]['start'],
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
        '{name}/{nbases}/merged/matrices/{name}.h5'
    log:
        'logs/homer2H5/{name}-{nbases}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    threads:
        workflow.cores
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 &> {log}'


rule mergeBins:
    input:
        rules.homer2H5.output
    output:
        '{name}/{nbases}/merged/matrices/{name}-{binsize}.h5'
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
            if config['HiC']['vMin2']:
                command += f' --vMin2 {config["HiC"]["vMin2"]} '
            if config['HiC']['vMax2']:
                command += f' --vMax2 {config["HiC"]["vMax2"]} '
        if config['genome']['genes']:
            command += f' --genes {config["genome"]["genes"]} '
        if config['ctcf']['data'] is not None:
            command += f'--ctcfOrient {rules.scaleCTCF.output} '
    if config['HiC']['vMin']:
        command += f' --vMin {config["HiC"]["vMin"]} '
    if config['HiC']['vMax']:
        command += f' --vMax {config["HiC"]["vMax"]} '
    if config['HiC']['log']:
        command += ' --logMatrix1 --logMatrix2 '

    command += f'--colourMap {config["HiC"]["colourMap"]}'
    return command


def getCTCFOrient(wc):
    """ Return CTCF orientation if not synthetic """
    if config['syntheticSequence'] is None:
        return rules.scaleCTCF.output
    else:
        return []

def getDepth(wc):
    return details[wc.name]['end'] - details[wc.name]['start'] + 1


rule createConfig:
    input:
        matrix = rules.mergeBins.output,
        ctcfOrient = getCTCFOrient
    output:
        '{name}/{nbases}/merged/config/{name}-{binsize}-configs.ini'
    params:
        depth = getDepth,
        hicConfig = getHiCconfig,
    group:
        'plotHiC'
    log:
        'logs/createConfig/{name}-{nbases}-{binsize}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input.matrix} '
        '--depth {params.depth} {params.hicConfig} > {output} 2> {log}'


def getTitle(wc):
    return f'"{wc.name} : {getRegion(wc)}"'


rule plotHiC:
    input:
        rules.createConfig.output
    output:
        'plots/contactMatrix/{name}-{nbases}-{binsize}-contactMatrix.png'
    params:
        region = getRegion,
        title = getTitle,
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
            '{name}/{nbases}/merged/matrices/contacts.pre.tsv'
        params:
            chr = lambda wc: details[wc.name]['chr'],
            start = lambda wc: details[wc.name]['start'],
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
            '{name}/{nbases}/merged/matrices/contacts.hic'
        log:
            'logs/juicerPre/{name}-{nbases}.log'
        params:
            chr = lambda wc: details[wc.name]['chr'],
            resolutions = '1000'
        resources:
             mem_mb = 16000
        threads:
            workflow.cores
        conda:
            f'{ENVS}/openjdk.yaml'
        shell:
            'java -Xmx{resources.mem_mb}m '
            '-jar {SCRIPTS}/juicer_tools_1.14.08.jar pre '
            '-c {params.chr} -r {params.resolutions} '
            '{input.tsv} {output} {input.chrom_sizes} &> {log}'


rule custom2XYZ:
    input:
        '{name}/{nbases}/reps/{rep}/lammps/{mode}.custom.gz'
    output:
        temp('{name}/{nbases}/reps/{rep}/lammps/{mode}.xyz.gz')
    group:
        'vmd'
    log:
        'logs/custom2XYZ/{name}-{nbases}-{rep}-{mode}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '(zcat {input} | {SCRIPTS}/custom2XYZ.py | gzip > {output}) 2> {log}'


checkpoint vmd:
    input:
        '{name}/{nbases}/reps/1/lammps/simulation.xyz.gz'
    output:
        directory('{name}/{nbases}/vmd/sequence')
    group:
        'vmd'
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


rule createGIF:
    input:
        rules.vmd.output,
        images = aggregateVMD
    output:
        'vmd/{name}-{nbases}-1-simulation.gif'
    params:
        delay = config['GIF']['delay'],
        loop = config['GIF']['loop']
    group:
        'vmd'
    log:
        'logs/createGIF/{name}-{nbases}.log'
    conda:
        f'{ENVS}/imagemagick.yaml'
    shell:
        'convert -delay {params.delay} -loop {params.loop} '
        '{input.images} {output} &> {log}'


## Archive ##

rule plotTUdistances:
    input:
        distances = expand(
            '{{name}}/{{nbases}}/reps/{rep}/TU-pairDistance.csv.gz', rep=REPS),
        beadDistribution = rules.writeTUdistribution.output
    output:
        '{name}/{nbases}/merged/{name}-TU-pairDistance.png',
    params:
        minRep = config['plotTU']['minRep'],
        fontSize = config['plotTU']['fontSize']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotTUdistances'
    log:
        'logs/plotTUdistances/{name}-{nbases}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotTUdistances.py {input.beadDistribution} '
        '{input.distances} --out {output} --fontSize {params.fontSize} '
        '--minRep {params.minRep} &> {log}'


def computeDTWparams(wc):
    """ Build paramters for computeDTW """

    sampletime = config[f'dtw{wc.type}']['sampletime']
    maxtime = config[f'dtw{wc.type}']['maxtime']
    params = f'--sampletime {sampletime} --maxtime {maxtime} '
    params += f'--group {wc.type} '

    return params


rule computeDTW:
    input:
        sim = rules.reformatLammps.output,
        groups = rules.getAtomGroups.output
    output:
        '{name}/{nbases}/reps/{rep}/{type}-euclideanDTW.csv.gz'
    params:
        computeDTWparams
    group:
        'computeAllDTW' if config['groupJobs'] else 'computeDTW'
    log:
        'logs/computeDTW/{name}-{nbases}-{rep}-{type}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/computeTUdtw.py --out {output} {params} '
        '{input.groups} {input.sim} &> {log}'


rule computeDTWnormalisation:
    input:
         expand('{{name}}/{{nbases}}/reps/{rep}/DNA-euclideanDTW.csv.gz',
            rep=REPS)
    output:
        '{name}/{nbases}/merged/info/DNA-normalisationFactors.csv.gz'
    log:
        'logs/computeDTWnormalisation/{name}-{nbases}.log'
    group:
        'computeAllDTW' if config['groupJobs'] else 'computeDTWnormalisation'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/computeDTWnormalisation.py --out {output} {input} &> {log}'


rule normaliseDTW:
    input:
        dtw = rules.computeDTW.output,
        normalisation = rules.computeDTWnormalisation.output
    output:
        '{name}/{nbases}/reps/{rep}/{type}-euclideanDTWnormalise.csv.gz'
    log:
        'logs/normaliseDTW/{name}-{nbases}-{rep}-{type}.log'
    group:
        'computeAllDTW' if config['groupJobs'] else 'normaliseDTW'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/normaliseDTW.py --out {output} '
        '{input.normalisation} {input.dtw} &> {log}'


rule plotDTW:
    input:
        dtw = expand(
            '{{name}}/{{nbases}}/reps/{rep}/{{type}}-euclideanDTWnormalise.csv.gz', rep=REPS),
        beadDistribution = rules.writeTUdistribution.output
    output:
        npz = '{name}/{nbases}/merged/matrices/{name}-{type}-dtwEuclidean.npz',
        plot = '{name}/{nbases}/merged/{name}-{type}-dtwEuclidean.png'
    params:
        minRep = config['plotTU']['minRep'],
        fontSize = config['plotTU']['fontSize']
    group:
        'computeAllDTW' if config['groupJobs'] else 'plotDTW'
    log:
        'logs/plotDTW/{name}-{nbases}-{type}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotDTWeuclidean.py {input.beadDistribution} '
        '{input.dtw} --out {output.plot} --npz {output.npz} '
        '--fontSize {params.fontSize} --minRep {params.minRep} &> {log}'


rule plotTUactivationByTime:
    input:
        rules.processTUinfo.output
    output:
        directory('{name}/{nbases}/plots/{rep}/')
    group:
        'plotTUactivationByTime'
    log:
        'logs/plotTUactivationByTime/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotTUactivationByTime.py --outdir {output} {input} &> {log}'


rule aggregateTarget:
    input:
        expand('{{name}}/{{nbases}}/plots/{rep}/', rep=REPS)
    output:
        touch(temp('{name}/{nbases}/plots/.aggregate.tmp'))
    group:
        'plotTUactivationByTime' if config['groupJobs'] else 'aggregateTarget'
