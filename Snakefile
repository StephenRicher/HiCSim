#!/usr/bin/env python3

container: "docker://continuumio/miniconda3:4.7.12"

import os
import math
import random
import tempfile
from itertools import combinations
from set_config import set_config, read_paths, getNbeads, adjustCoordinates

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
    'ctcf':           None,
    'masking':        {},
    'genome' :        {'build':     'genome',
                       'chromSizes': None,
                       'name':       None,
                       'genes':      None,
                       'chr':        None,
                       'start':      None,
                       'end':        None,},
    'retentionProb':   [0.5],
    'syntheticSequence' : {}              ,
    'basesPerBead':   2000,
    'monomers':       100,
    'method':         'mean',
    'coeffs':         '',
    'reps':           5,
    'maxReps':        0,
    'equilibrateOnly': False,
    'extraStats':      True,
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
    'lammps':         {'restart':          0    ,
                       'writeInterval':    10   ,
                       'timestep':         0.01 ,
                       'warmUp':           20000,
                       'simTime':          20000,
                       'harmonicCoeff':    2    ,
                       'noExtrusion':      False,
                       'TFswap':           100  ,
                       'extrudersPerMb':   8,
                       'nSplit':           10   ,
                       'threads':          1    ,},
    'HiC':            {'matrix' :    None    ,
                       'chrom':      None    ,
                       'log' :       True    ,
                       'colourMap': 'Purples',
                       'dpi':        300     ,
                       'vMin':       None    ,
                       'vMax':       None    ,
                       'vMin2':      None    ,
                       'vMax2':      None    ,},
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
}
config = set_config(config, default_config)

workdir : config['workdir']

details = {}
if not config['syntheticSequence']:
    # If no synthetic sequence ensure genome information is provided.
    invalid = False
    for key in ['chromSizes', 'name', 'chr', 'start', 'end']:
        if config['genome'][key] is None:
            sys.stderr.write(
                f'\033[31mNo synthetic sequence provided and '
                f'config["genome"]["{key}"] not set.\033[m\n')
            invalid = True
    if invalid:
        sys.exit('\033[31mInvalid configuration setting.\033[m\n')
    for name in [config['genome']['name']]:
        start, end, nBeads = adjustCoordinates(
            config['genome']['start'],
            config['genome']['end'],
            config['basesPerBead'])
        scaledRegion = f"{config['genome']['chr']}-{start}-{end}"
        print(f'Adjusting {name} positions to {scaledRegion}', file=sys.stderr)
        details[name] = {'chr':   config['genome']['chr'],
                         'start': start,
                         'end': end,
                         'nBeads': nBeads}
    # Subtraction matrix not yet supported for 'real' sequences
    COMPARES = []
else:
    config['genome']['chromSizes'] = []
    allBeadLengths = []
    for name in config['syntheticSequence'].keys():
        nBeads = getNbeads(config['syntheticSequence'][name])
        allBeadLengths.append(nBeads)
        end = (nBeads * config['basesPerBead'])
        details[name] = {'chr':    'synthetic',
                         'start':  1,
                         'end':    end,
                         'nBeads': nBeads}
    if len(set(allBeadLengths)) != 1:
        print('Subtraction matrices only supported if all synthetic '
              'sequences are the same length', file=sys.stderr)
        COMPARES = []
    else:
        COMPARES = [f'{i[0]}-vs-{i[1]}' for i in combinations(details, 2)]

BUILD = config['genome']['build']

# Define list of N reps from 1 to N
REPS = list(range(1, config['reps'] + 1))

MAXREPS = list(range(1, max(config['maxReps'], config['reps']) + 1))

# Set seeds for different parts of workflow.
seeds = {}
for type in ['sequence', 'initialConform', 'monomerPositions', 'simulation']:
    random.seed(config['random']['seed'])
    if config['random'][type]:
        seeds[type] = [random.randint(1, (2**16) - 1) for rep in MAXREPS]
    else:
        seeds[type] = [random.randint(1, (2**16) - 1)] * len(MAXREPS)

if config['lammps']['simTime'] < config['lammps']['writeInterval']:
    sys.exit('Simulation time less than writeInterval')

if (workflow.cores is not None) and (config['lammps']['threads'] > workflow.cores):
    print(f'Lammps threads cannot be higher than provided workflow cores. '
          f'Adjusting lammps threads from {config["lammps"]["threads"]} '
          f'to {workflow.cores}.', file=sys.stderr)
    config['lammps']['threads'] = workflow.cores

# Read track file basename and masking character into a dictionary
track_data = {}
if config['masking']:
    for file, character in config['masking'].items():
        track_file = f'{os.path.basename(file)}'
        track_data[track_file] = {'source' : file, 'character' : character}

wildcard_constraints:
    all = r'[^\/]+',
    stat = r'stats|pairStats|TADstatus',
    name = rf'{"|".join(details.keys())}',
    nbases = rf'{config["basesPerBead"]}',
    rep = rf'{"|".join([str(rep) for rep in MAXREPS])}'

extraStats = ([
    expand('plots/{plot}/{name}-{nbases}-{prob}-{rep}-{plot}.png',
        nbases=config['basesPerBead'], name=details.keys(), rep=REPS,
        plot=['TADstructure'], prob=config['retentionProb']),
    expand('plots/{plot}/{name}-{nbases}-{prob}-{plot}.png',
        nbases=config['basesPerBead'], name=details.keys(),
        plot=['meanVariance','TUcorrelation', 'TUcircos',
              'radiusGyration', 'TUreplicateCount'],
        prob=config['retentionProb']),
    expand('plots/pairCluster/{name}-{nbases}-{prob}.png',
        nbases=config['basesPerBead'], name=details.keys(),
        prob=config['retentionProb']),
    expand('{name}/{nbases}/merged/{name}-TU-{stat}-{prob}.csv.gz',
        name=details.keys(), nbases=config['basesPerBead'],
        stat=['stats', 'pairStats', 'TADstatus'],
        prob=config['retentionProb']),
]) if config['extraStats'] else []

rule all:
    input:
        extraStats,
        expand('{name}/{nbases}/lammpsInit/simulation-equil-{prob}',
            name=details.keys(), nbases=config['basesPerBead'],
            prob=config['retentionProb']),
        ([expand('vmd/{name}-{nbases}-{prob}-1-simulation.gif',
            nbases=config['basesPerBead'], prob=config['retentionProb'],
            name=details.keys()) if config['GIF']['create'] else [],
         expand('plots/contactMatrix/{name}-{nbases}-{prob}-contactMatrix.svg',
            nbases=config['basesPerBead'], name=details.keys(),
            prob=config['retentionProb']),
         expand('plots/radiusGyration/{name}-{nbases}-{prob}-radiusGyration.png',
            nbases=config['basesPerBead'], name=details.keys(),
            prob=config['retentionProb']),
        expand('{name}/{nbases}/reps/{rep}/lammps/config/atomGroups-{monomers}-{prob}.json',
            name=details.keys(), nbases=config['basesPerBead'],
            rep=REPS, monomers=config['monomers'],
            prob=config['retentionProb']),
        expand('comparison/{nbases}/{compare}-{nbases}-{prob}-logFC.svg',
            nbases=config['basesPerBead'], compare=COMPARES,
            prob=config['retentionProb']),
        (expand('HiCRep/{name}/{name}-{nbases}-{prob}-vs-experimental.csv',
            name=details.keys(), nbases=config['basesPerBead'],
            prob=config['retentionProb'])
            if config['HiC']['matrix'] else [])]
        if not config['equilibrateOnly'] else [])


if config['ctcf'] is not None:

    rule filterCTCF:
        input:
            config['ctcf']
        output:
            'tracks/{rep}/CTCF/CTCF-filtered-{prob}.bed.gz'
        params:
            seed = lambda wc: seeds['sequence'][int(wc.rep) - 1]
        group:
            'prepLammps'
        log:
            'logs/sampleBed/{rep}-{prob}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            '({SCRIPTS}/filterBedScore.py {input} F R --seed {params.seed} '
            '--retentionProb {wildcards.prob} | gzip > {output}) 2> {log}'


rule filterTracks:
    input:
        lambda wc: track_data[wc.track]['source']
    output:
        'tracks/{rep}/other/{track}-{prob}.gz'
    params:
        seed = lambda wc: seeds['sequence'][int(wc.rep) - 1],
        character = lambda wc: track_data[wc.track]['character']
    group:
        'prepLammps'
    log:
        'logs/filterTracks/{track}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '({SCRIPTS}/filterBedScore.py {input} {params.character} '
        '--seed {params.seed} --retentionProb {wildcards.prob} '
        '| gzip > {output}) 2> {log}'


def getRegion(wc):
    try:
        name = wc.name
    except AttributeError:
        name = wc.name1
    chrom = details[name]['chr']
    start = details[name]['start']
    end = details[name]['end']
    return f'{chrom}:{start}-{end}'


def getCTCF(wc):
    """ Return rule output only if used. """
    if config['ctcf'] is not None:
        return rules.filterCTCF.output
    else:
        return []


rule maskFasta:
    input:
        getCTCF,
        expand('tracks/{{rep}}/other/{track}-{{prob}}.gz', track=track_data.keys()),
    output:
        '{name}/{nbases}/reps/{rep}/{name}-{prob}-beads-{rep}.txt'
    params:
        region = getRegion,
        nBases = config['basesPerBead'],
        seed = lambda wc: seeds['sequence'][int(wc.rep) - 1]
    group:
        'prepLammps'
    log:
        'logs/maskFasta/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/maskFasta.py {params.region} {input} '
        '--seed {params.seed} --nBases {params.nBases} '
        '> {output} 2> {log}'


rule sampleSynthetic:
    input:
        lambda wc: config['syntheticSequence'][wc.name]
    output:
        '{name}/{nbases}/reps/{rep}/sampledSynthetic-{rep}.txt'
    group:
        'lammpsEquilibrate'
    params:
        seed = lambda wc: seeds['sequence'][int(wc.rep) - 1],
    log:
        'logs/sampleSynthetic/{name}-{nbases}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/sampleSynthetic.py {input} --seed {params.seed} '
        '> {output} 2> {log}'


def allBeadsInput(wc):
    if config['syntheticSequence']:
        return expand(
            '{{name}}/{{nbases}}/reps/{rep}/sampledSynthetic-{rep}.txt', rep=MAXREPS)
    else:
        return expand('{{name}}/{{nbases}}/reps/{rep}/{{name}}-{{prob}}-beads-{rep}.txt', rep=MAXREPS)


rule extractAtomTypes:
    input:
        allBeadsInput
    output:
        '{name}/{nbases}/beadTypeID-{prob}.json'
    group:
        'lammpsEquilibrate'
    log:
        'logs/extractAtomTypes/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/extractAtomTypes.py {input} > {output} 2> {log}'


def setBeadsToLammpsCmd():
    cmd = '{SCRIPTS}/generate_polymer.py '
    cmd += (
        '--polymerSeed {params.polymerSeed} '
        '--monomerSeed {params.monomerSeed} '
        '--xlo {params.xlo} --xhi {params.xhi} '
        '--ylo {params.ylo} --yhi {params.yhi} '
        '--zlo {params.zlo} --zhi {params.zhi} '
        '--nMonomers {params.nMonomers} '
        '--basesPerBead {params.basesPerBead} '
        '--beadTypes {input} '
        '{params.randomWalk} {params.nBeads} > {output} 2> {log}')
    return cmd


rule BeadsToLammps:
    input:
        rules.extractAtomTypes.output
    output:
        f'{{name}}/{{nbases}}/lammpsInit/lammps_input-{config["monomers"]}-{{prob}}.dat'
    params:
        nMonomers = config['monomers'],
        nBeads = lambda wc: details[wc.name]['nBeads'],
        basesPerBead = config['basesPerBead'],
        polymerSeed = 2, #lambda wc: seeds['initialConform'][int(wc.rep) - 1],
        monomerSeed = 2, #lambda wc: seeds['monomerPositions'][int(wc.rep) - 1],
        xlo = config['box']['xlo'],
        xhi = config['box']['xhi'],
        ylo = config['box']['ylo'],
        yhi = config['box']['yhi'],
        zlo = config['box']['zlo'],
        zhi = config['box']['zhi'],
        randomWalk = '--randomWalk' if config['random']['walk'] else ''
    group:
        'lammpsEquilibrate'
    log:
        'logs/BeadsToLammps/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        setBeadsToLammpsCmd()

def setLmpPrefix(wc):
    """ Define Lammps command for mpi or not """
    if config['lammps']['threads'] > 1:
        return f'mpirun -np {config["lammps"]["threads"]}'
    else:
        return ''


rule lammpsEquilibrate:
    input:
        rules.BeadsToLammps.output
    output:
        equil = '{name}/{nbases}/lammpsInit/simulation-equil-{prob}',
        equilInfo = '{name}/{nbases}/lammpsInit/warmUp-{prob}.custom.gz',
        radiusGyration = '{name}/{nbases}/lammpsInit/radiusOfGyration-{prob}.txt'
    params:
        writeInterval = config['lammps']['writeInterval'],
        timestep = config['lammps']['timestep'],
        seed = 1,
        cosinePotential = lambda wc: 10000 / config['basesPerBead'],
        equilTime = config['lammps']['warmUp'],
        lmpPrefix = setLmpPrefix
    group:
        'lammpsEquilibrate'
    log:
        'logs/lammpsEquilibrate/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    threads:
        config['lammps']['threads']
    shell:
        '{params.lmpPrefix} python {SCRIPTS}/runLammpsEquilibrate.py {input} '
        '{output.equil} --timestep {params.timestep} '
        '--writeInterval {params.writeInterval} --equilTime {params.equilTime} '
        '--seed {params.seed} --equilInfo {output.equilInfo} '
        '--cosinePotential {params.cosinePotential} '
        '--radiusGyrationOut {output.radiusGyration} &> {log}'


def beadsInput(wc):
    if config['syntheticSequence']:
        return rules.sampleSynthetic.output
    else:
        return rules.maskFasta.output


rule getAtomGroups:
    input:
        beadsInput
    output:
        f'{{name}}/{{nbases}}/reps/{{rep}}/lammps/config/atomGroups-{config["monomers"]}-{{prob}}.json'
    params:
        TUs = ['P','p'],
        nMonomers = config['monomers'],
    group:
        'processAllLammps' if config['groupJobs'] else 'getAtomGroups'
    log:
        'logs/getAtomGroups/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/getAtomGroups.py {input} --TUs {params.TUs} '
        '--nMonomers {params.nMonomers} > {output} 2> {log}'


rule lammpsSimulation:
    input:
        equil = rules.lammpsEquilibrate.output.equil,
        groups = rules.getAtomGroups.output,
        beadTypes = rules.extractAtomTypes.output
    output:
        simOut = '{name}/{nbases}/reps/{rep}/lammps/simulation-{prob}.custom.gz',
        TADstatus = '{name}/{nbases}/reps/{rep}/beadTADstatus-{prob}.csv.gz',
        radiusGyration = '{name}/{nbases}/reps/{rep}/lammps/radiusOfGyration-{prob}.txt',
    params:
        simTime = config['lammps']['simTime'],
        noExtrusion = '--noExtrusion' if config['lammps']['noExtrusion'] else '',
        TFswap = config['lammps']['TFswap'],
        coeffs = config['coeffs'],
        timestep = config['lammps']['timestep'],
        nBasesPerBead = config['basesPerBead'],
        extrudersPerMb = config['lammps']['extrudersPerMb'],
        writeInterval = config['lammps']['writeInterval'],
        seed = lambda wc: seeds['simulation'][int(wc.rep) - 1],
        harmonicCoeff = config['lammps']['harmonicCoeff'],
        sim = '{name}/{nbases}/reps/{rep}/lammps/simulation.custom.gz',
        lmpPrefix = setLmpPrefix
    group:
        'processAllLammps' if config['groupJobs'] else 'lammpsSimulation'
    log:
        'logs/lammps/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/lammps.yaml'
    threads:
        config['lammps']['threads']
    shell:
        '({params.lmpPrefix} python {SCRIPTS}/runLammpsSimulation.py --verbose '
        '{input.equil} {input.groups} '
        '--simTime {params.simTime} --writeInterval {params.writeInterval} '
        '--seed {params.seed} --harmonicCoeff {params.harmonicCoeff} '
        '--TFswap {params.TFswap} --radiusGyrationOut {output.radiusGyration} '
        '--simOut {output.simOut} --timestep {params.timestep} '
        '{params.noExtrusion} --nBasesPerBead {params.nBasesPerBead} '
        '--pairCoeffs {params.coeffs} --beadTypes {input.beadTypes} '
        '--nExtrudersPerMb {params.extrudersPerMb} | gzip > {output.TADstatus}) &> {log}'


rule plotRG:
    input:
        '{name}/{nbases}/lammpsInit/radiusOfGyration-{prob}.txt',
        expand('{{name}}/{{nbases}}/reps/{rep}/lammps/radiusOfGyration-{{prob}}.txt',
            rep=REPS)
    output:
        'plots/radiusGyration/{name}-{nbases}-{prob}-radiusGyration.png'
    params:
        confidence = config['plotRG']['confidence'],
        dpi = config['plotRG']['dpi']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotRG'
    log:
        'logs/plotRG/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotRG.py --out {output} --dpi {params.dpi} '
        '--confidence {params.confidence} {input} 2> {log}'


rule plotTADstucture:
    input:
        rules.lammpsSimulation.output.TADstatus
    output:
        'plots/TADstructure/{name}-{nbases}-{prob}-{rep}-TADstructure.png'
    params:
        dpi = config['plotRG']['dpi']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotTADstucture'
    log:
        'logs/plotTADstucture/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/plotTADstructure.py --dpi {params.dpi} '
        '--out {output} {input} &> {log}'


rule writeTUdistribution:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/lammps/config/atomGroups-{monomers}-{{prob}}.json',
            rep=REPS, monomers=config["monomers"])
    output:
        '{name}/{nbases}/merged/info/{name}-{prob}-TU-distribution.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'writeTUdistribution'
    log:
        'logs/writeTUdistribution/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/writeTUdistribution.py --out {output} {input} &> {log}'


def getNsteps(wc):
    """ Get number of simulation steps to write
        per file to generate nSplit files"""
    # Compute number of simulation write times
    nWrites = math.floor(config['lammps']['simTime'] / config['lammps']['writeInterval']) + 1
    # Number of steps per file to generate nSplit files
    nSteps = math.ceil(nWrites / config['lammps']['nSplit'])
    return nSteps


rule reformatLammps:
    input:
        rules.lammpsSimulation.output.simOut
    output:
        temp(expand('{{name}}/{{nbases}}/reps/{{rep}}/lammps/simulation-{{prob}}-split{split}.csv.gz',
            split=range(config['lammps']['nSplit'])))
    params:
        nSteps = getNsteps,
        prefix = lambda wc: f'{wc.name}/{wc.nbases}/reps/{wc.rep}/lammps/simulation-{wc.prob}-split'
    group:
        'processAllLammps' if config['groupJobs'] else 'reformatLammps'
    log:
        'logs/reformatLammps/{name}-{nbases}-{prob}-{rep}.log'
    shell:
        '{SCRIPTS}/reformatLammps.awk -v nSteps={params.nSteps} '
        '-v prefix={params.prefix} <(zcat {input}) &> {log}'


rule processTUinfo:
    input:
        atomGroups = rules.getAtomGroups.output,
        TADstatus = rules.lammpsSimulation.output.TADstatus,
        sim = '{name}/{nbases}/reps/{rep}/lammps/simulation-{prob}-split{split}.csv.gz',
    output:
        '{name}/{nbases}/reps/{rep}/TU-info-{prob}-split{split}.csv.gz',
    params:
        distance = 1.8
    group:
        'processAllLammps' if config['groupJobs'] else 'processTUinfo'
    log:
        'logs/processTUinfo/{name}-{nbases}-{prob}-{rep}-{split}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/processTUinfo.py --out {output} '
        '--distance {params.distance} {input.atomGroups} '
        '{input.TADstatus} {input.sim} &> {log}'


rule mergeTUinfo:
    input:
        expand('{{name}}/{{nbases}}/reps/{{rep}}/TU-info-{{prob}}-split{split}.csv.gz',
            split=range(config['lammps']['nSplit'])),
    output:
        '{name}/{nbases}/reps/{rep}/TU-info-{prob}.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'processTUinfo'
    log:
        'logs/mergeTUinfo/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mergeByRep.py --noRep --out {output} {input} &> {log}'


rule processTADstatus:
    input:
        rules.mergeTUinfo.output
    output:
        '{name}/{nbases}/reps/{rep}/TU-TADstatus-{prob}.csv.gz',
    group:
        'processAllLammps' if config['groupJobs'] else 'processTUinfo'
    log:
        'logs/processTADstatus/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/processTADstatus.py {input} --out {output} &> {log}'


rule DBSCAN:
    input:
        rules.mergeTUinfo.output
    output:
        clusterPairs = '{name}/{nbases}/reps/{rep}/TU-clusterPair-{prob}.csv.gz',
        clusterPlot = 'plots/DBSCAN/{name}/{name}-{nbases}-{prob}-{rep}-cluster.png',
    params:
        eps = 6,
        minSamples = 2
    group:
        'DBSCAN_all' if config['groupJobs'] else 'DBSCAN'
    log:
        'logs/DBSCAN/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/dbscanTU.py --outplot {output.clusterPlot} '
        '--eps {params.eps} --minSamples {params.minSamples} '
        '--out {output.clusterPairs} {input} &> {log}'


rule plotDBSCAN:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/TU-clusterPair-{{prob}}.csv.gz', rep=REPS),
    output:
         'plots/pairCluster/{name}-{nbases}-{prob}-pairCluster.png'
    group:
        'DBSCAN_all' if config['groupJobs'] else 'plotDBSCAN'
    log:
        'logs/plotDBSCAN/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotDBSCAN.py {output} {input} &> {log}'


rule computeTUstats:
    input:
        rules.mergeTUinfo.output
    output:
        TUstats = '{name}/{nbases}/reps/{rep}/TU-stats-{prob}.csv.gz',
        TUpairStats = '{name}/{nbases}/reps/{rep}/TU-pairStats-{prob}.csv.gz',
    group:
        'processAllLammps' if config['groupJobs'] else 'computeTUstats'
    log:
        'logs/computeTUstats/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/computeTUstats.py --out {output.TUstats} '
        '--TUpairStats {output.TUpairStats} {input} &> {log}'


rule mergeByRep:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/TU-{{stat}}-{{prob}}.csv.gz', rep=REPS),
    output:
        '{name}/{nbases}/merged/{name}-TU-{stat}-{prob}.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'mergeByRep'
    log:
        'logs/mergeByRep/{name}-{nbases}-{prob}-{stat}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mergeByRep.py --out {output} {input} &> {log}'


rule plotMeanVariance:
    input:
        '{name}/{nbases}/merged/{name}-TU-stats-{prob}.csv.gz'
    output:
        'plots/meanVariance/{name}-{nbases}-{prob}-meanVariance.png'
    params:
        fontSize = config['plotTU']['fontSize']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotMeanVariance'
    log:
        'logs/plotMeanVariance/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotMeanVariance.py {input} --out {output} '
        '--fontSize {params.fontSize}  &> {log}'


rule computeTUcorrelation:
    input:
        '{name}/{nbases}/merged/{name}-TU-stats-{prob}.csv.gz'
    output:
        '{name}/{nbases}/merged/{name}-TU-correlation-{prob}.csv.gz'
    group:
        'processAllLammps' if config['groupJobs'] else 'computeTUcorrelation'
    log:
        'logs/computeTUcorrelation/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/computeTUcorrelation.py --out {output} {input} &> {log}'


rule plotTUcorrelation:
    input:
        correlations = rules.computeTUcorrelation.output,
        beadDistribution = rules.writeTUdistribution.output
    output:
        meanHeatmap = 'plots/TUcorrelation/{name}-{nbases}-{prob}-TUcorrelation.png',
        sumHeatmap = 'plots/TUreplicateCount/{name}-{nbases}-{prob}-TUreplicateCount.png',
        circos = 'plots/TUcircos/{name}-{nbases}-{prob}-TUcircos.png'
    params:
        pvalue = config['plotTU']['pvalue'],
        vMin = config['plotTU']['vMin'],
        vMax = config['plotTU']['vMax'],
        minRep = config['plotTU']['minRep'],
        fontSize = config['plotTU']['fontSize']
    group:
        'processAllLammps' if config['groupJobs'] else 'plotTUcorrelation'
    log:
        'logs/plotTUcorrelation/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/plotTUcorrelation.py {input.beadDistribution} '
        '{input.correlations} --circos {output.circos} '
        '--sumHeatmap {output.sumHeatmap} --meanHeatmap {output.meanHeatmap} '
        '--fontSize {params.fontSize} --pvalue {params.pvalue} '
        '--minRep {params.minRep} --vmin {params.vMin} '
        '--vmax {params.vMax} &> {log}'


rule createContactMatrix:
    input:
        xyz = '{name}/{nbases}/reps/{rep}/lammps/simulation-{prob}-split{split}.csv.gz',
        groups = rules.getAtomGroups.output
    output:
        '{name}/{nbases}/reps/{rep}/matrices/contacts-{prob}-split{split}.npz'
    params:
        distance = 9 / (config['basesPerBead'] / 1000),
        periodic = '--periodic',
        seed = lambda wc: seeds['simulation'][int(wc.rep) - 1],
        x = abs(config['box']['xhi'] - config['box']['xlo']),
        y = abs(config['box']['yhi'] - config['box']['ylo']),
        z = abs(config['box']['zhi'] - config['box']['zlo'])
    group:
        'processAllLammps' if config['groupJobs'] else 'createContactMatrix'
    log:
        'logs/createContactMatrix/{name}-{nbases}-{prob}-{rep}-{split}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/createContactMatrix.py {params.periodic} '
        '--out {output} --distance {params.distance} --seed {params.seed} '
        '--dimensions {params.x} {params.y} {params.z} {input.groups} '
        '<(zcat -f {input.xyz}) &> {log}'


rule mergeSplit:
    input:
        expand('{{name}}/{{nbases}}/reps/{{rep}}/matrices/contacts-{{prob}}-split{split}.npz',
            split=range(config['lammps']['nSplit']))
    output:
        '{name}/{nbases}/reps/{rep}/matrices/contacts-{prob}.npz'
    params:
        method = 'sum'
    group:
        'processAllLammps' if config['groupJobs'] else 'createContactMatrix'
    log:
        'logs/mergeSplit/{name}-{nbases}-{prob}-{rep}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mergeMatrices.py --out {output} '
        '--method {params.method} {input} &> {log}'


rule mergeReplicates:
    input:
        expand('{{name}}/{{nbases}}/reps/{rep}/matrices/contacts-{{prob}}.npz',
            rep=REPS)
    output:
        '{name}/{nbases}/merged/matrices/{name}-{prob}.npz'
    params:
        method = config['method']
    group:
        'processAllLammps' if config['groupJobs'] else 'mergeReplicates'
    log:
        'logs/mergeReplicates/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/mergeMatrices.py --out {output} '
        '--method {params.method} {input} &> {log}'


def set_expHiC(wc):
    if config['HiC']['matrix']:
        return f'--expHiC {config["HiC"]["matrix"]}'
    else:
        return ''


rule matrix2homer:
    input:
        '{name}/{nbases}/merged/matrices/{name}-{prob}.npz'
    output:
        '{name}/{nbases}/merged/matrices/{name}-{prob}.homer'
    params:
        chr = lambda wc: details[wc.name]['chr'],
        start = lambda wc: details[wc.name]['start'],
        binSize = config['basesPerBead'],
        expHiC = set_expHiC
    log:
        'logs/matrix2homer/{name}-{nbases}-{prob}log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/npz2homer.py --chromosome {params.chr} '
        '{params.expHiC} --start {params.start} '
        '--binSize {params.binSize} {input} > {output} 2> {log}'


rule homer2H5:
    input:
        rules.matrix2homer.output
    output:
        '{name}/{nbases}/merged/matrices/{name}-{prob}.h5'
    log:
        'logs/homer2H5/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicConvertFormat --matrices {input} --outFileName {output} '
        '--inputFormat homer --outputFormat h5 &> {log}'


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
        if config['ctcf'] is not None:
            command += f'--ctcfOrient {config["ctcf"]} '
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
        return config["ctcf"]
    else:
        return []

def getDepth(wc):
    try:
        name = wc.name
    except AttributeError:
        name = wc.name1
    return details[name]['end'] - details[name]['start'] + 1


rule createConfig:
    input:
        matrix = rules.homer2H5.output,
        ctcfOrient = getCTCFOrient,
    output:
        '{name}/{nbases}/merged/config/{name}-{prob}-configs.ini'
    params:
        depth = getDepth,
        hicConfig = getHiCconfig,
    group:
        'plotHiC'
    log:
        'logs/createConfig/{name}-{nbases}-{prob}.log'
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
        'plots/contactMatrix/{name}-{nbases}-{prob}-contactMatrix.svg'
    params:
        region = getRegion,
        title = getTitle,
        dpi = config['HiC']['dpi']
    group:
        'plotHiC'
    log:
        'logs/plotHiC/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    shell:
        'pyGenomeTracks --tracks {input} --region {params.region} '
        '--outFileName {output} --title {params.title} --dpi {params.dpi} '
        '&> {log}'


rule distanceNormalise:
    input:
        '{name}/{nbases}/merged/matrices/{name}-{prob}.h5'
    output:
        '{name}/{nbases}/merged/matrices/{name}-{prob}-obsExp.h5'
    log:
        'logs/distanceNormalise/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'hicTransform -m {input} --method obs_exp -o {output} &> {log}'


rule compareMatrices:
    input:
        '{name1}/{nbases}/merged/matrices/{name1}-{prob}.h5',
        '{name2}/{nbases}/merged/matrices/{name2}-{prob}.h5'
    output:
        'comparison/{nbases}/{name1}-vs-{name2}-{nbases}-{prob}.h5'
    log:
        'logs/compareMatrices/{name1}-vs-{name2}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/hicexplorer.yaml'
    shell:
        'python {SCRIPTS}/compareHiC.py {input} --outMatrix {output} &> {log}'


rule createCompareConfig:
    input:
        rules.compareMatrices.output
    output:
        'comparison/{nbases}/{name1}-vs-{name2}-{nbases}-{prob}-logFC.ini'
    params:
        depth = getDepth,
    log:
        'logs/creatCompareConfig/{name1}-vs-{name2}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '{SCRIPTS}/generate_config.py --matrix {input} --depth {params.depth} '
        '--vMin -2 --vMax 2 --colourMap bwr > {output} 2> {log}'


def getCompareTitle(wc):
    return f'"{wc.name1} vs. {wc.name2}: {getRegion(wc)}"'


rule plotSubtractionMatrix:
    input:
        rules.createCompareConfig.output
    output:
        'comparison/{nbases}/{name1}-vs-{name2}-{nbases}-{prob}-logFC.svg'
    params:
        region = getRegion,
        title = getCompareTitle,
        dpi = config['HiC']['dpi']
    log:
        'logs/plotSubtractionMatrix/{name1}-vs-{name2}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/pygenometracks.yaml'
    shell:
        'pyGenomeTracks --tracks {input} --region {params.region} '
        '--outFileName {output} --title {params.title} --dpi {params.dpi} '
        '&> {log}'


rule H5_to_NxN3p:
    input:
        rules.homer2H5.output
    output:
        '{name}/{nbases}/merged/matrices/{name}-{prob}.nxn3p.tsv'
    log:
        'logs/H5_to_NxN3p/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        'python {SCRIPTS}/H5_to_NxN3p.py {input} > {output} 2> {log}'


if config['HiC']['matrix'] is not None:

    rule H5_to_NxN3p_exp:
        input:
            config['HiC']['matrix']
        output:
            'expMatrix/{name}.nxn3p.tsv'
        params:
            chrom = lambda wc: details[wc.name]['chr'],
            start = lambda wc: details[wc.name]['start'],
            end = lambda wc: details[wc.name]['end'],
        log:
            'logs/H5_to_NxN3p/{name}.log'
        conda:
            f'{ENVS}/python3.yaml'
        shell:
            'python {SCRIPTS}/H5_to_NxN3p.py --chrom {params.chrom} '
            '--start {params.start} --end {params.end} {input} '
            '> {output} 2> {log}'


    rule HiCRep:
        input:
            'expMatrix/{name}.nxn3p.tsv',
            '{name}/{nbases}/merged/matrices/{name}-{prob}.nxn3p.tsv',
        output:
            'HiCRep/{name}/{name}-{nbases}-{prob}-vs-experimental.csv'
        params:
            start = lambda wc: details[wc.name]['start'],
            end = lambda wc: details[wc.name]['end']
        log:
            'logs/HiCRep/{name}-{nbases}-{prob}-vs-experimental.log'
        conda:
            f'{ENVS}/hicrep.yaml'
        shell:
            'Rscript {SCRIPTS}/runHiCRep.R 10000 '
            '{params.start} {params.end} {input} > {output} 2> {log}'


if config['syntheticSequence'] is None:

    rule matrix2pre:
        input:
            rules.mergeReplicates.output
        output:
            '{name}/{nbases}/merged/matrices/contacts.pre.tsv'
        params:
            chr = lambda wc: details[wc.name]['chr'],
            start = lambda wc: details[wc.name]['start'],
            binsize = config['basesPerBead']
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
            chrom_sizes = config['genome']['chromSizes']
        output:
            '{name}/{nbases}/merged/matrices/contacts.hic'
        log:
            'logs/juicerPre/{name}-{nbases}.log'
        params:
            chr = lambda wc: details[wc.name]['chr'],
            resolutions = '1000'
        resources:
             mem_mb = 16000
        conda:
            f'{ENVS}/openjdk.yaml'
        shell:
            'java -Xmx{resources.mem_mb}m '
            '-jar {SCRIPTS}/juicer_tools_1.14.08.jar pre '
            '-c {params.chr} -r {params.resolutions} '
            '{input.tsv} {output} {input.chrom_sizes} &> {log}'


rule custom2XYZ:
    input:
        '{name}/{nbases}/reps/{rep}/lammps/{mode}-{prob}.custom.gz'
    output:
        temp('{name}/{nbases}/reps/{rep}/lammps/{mode}-{prob}.xyz.gz')
    group:
        'vmd'
    log:
        'logs/custom2XYZ/{name}-{nbases}-{prob}-{rep}-{mode}.log'
    conda:
        f'{ENVS}/python3.yaml'
    shell:
        '(zcat {input} | {SCRIPTS}/custom2XYZ.py | gzip > {output}) 2> {log}'


checkpoint vmd:
    input:
        '{name}/{nbases}/reps/1/lammps/simulation-{prob}.xyz.gz'
    output:
        directory('{name}/{nbases}/vmd/sequence-{prob}')
    group:
        'vmd'
    log:
        'logs/vmd/{name}-{nbases}-{prob}.log'
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
        'vmd/{name}-{nbases}-{prob}-1-simulation.gif'
    params:
        delay = config['GIF']['delay'],
        loop = config['GIF']['loop']
    group:
        'vmd'
    log:
        'logs/createGIF/{name}-{nbases}-{prob}.log'
    conda:
        f'{ENVS}/imagemagick.yaml'
    shell:
        'convert -delay {params.delay} -loop {params.loop} '
        '{input.images} {output} &> {log}'
