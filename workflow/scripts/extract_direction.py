#!/usr/bin/env python3

import re
import sys
import subprocess
import pyCommonTools as pct


genome = sys.argv[1]
file = sys.argv[2]
chr_ref = sys.argv[3]
#start_ref = int(sys.argv[4])
#end_ref = int(sys.argv[5])
start_ref = None
end_ref = None

ctcf_seqs = {'+' : '>ctcf\nCCACNAGGTGGCAG',
            '-' : '>ctcf\nCTGCCACCTNGTGG'}


def within_region(chr, start, end, chr_ref=None, start_ref=None, end_ref=None):

    log = pct.create_logger()
    if chr_ref == None:
        rc = True
    elif chr_ref == chr:
        if start_ref == end_ref == None:
            rc = True
        elif start_ref != end_ref and (start_ref == None or end_ref == None):
            log.error('Only 1 of start and end coordinates are set.')
            sys.exit(1)
        elif start_ref > end_ref:
            log.error(
                f'Start position {start} is larger than end position {end}.')
            sys.exit(1)
        elif start >= start_ref and end <= end_ref:
            rc = True
        else:
            rc = False
    else:
        rc = False

    return rc


with open(file) as f:
    for record in f:
        if record.startswith('#'):
            continue

        record = record.strip().split()
        chr = record[1].strip('chr')
        start = int(record[2])
        end = int(record[3])

        if not within_region(chr, start, end, chr_ref, start_ref, end_ref):
            continue

        coordinates = f'{chr}:{start}-{end}'

        faidx_command = ['samtools', 'faidx', genome, coordinates]
        faidx = subprocess.Popen(faidx_command,
            stdout = subprocess.PIPE, stderr = subprocess.DEVNULL)
        genomic_seq = faidx.communicate()[0]

        if faidx.returncode != 0:
            continue

        alignments = {'+' : {}, '-' : {}}

        for orientation in ctcf_seqs:
            ctcf_seq = ctcf_seqs[orientation]

            needle = subprocess.Popen(
                [f'{sys.path[0]}/run_needle.sh', genomic_seq, ctcf_seq],
                stdout = subprocess.PIPE, stderr = subprocess.DEVNULL)

            ctcf_alignment_line = ''
            for line in needle.stdout:
                line = line.decode('ascii')
                if line.startswith('# Score'):
                    alignments[orientation]['score'] = float(line.split()[-1])
                elif line.startswith('ctcf'):
                    ctcf_alignment_line += line.split()[2]


            ctcf_min = ctcf_max = None
            for index, char in enumerate(ctcf_alignment_line):
                if char != '-':
                    if ctcf_min is None:
                        ctcf_min = index
                    if ctcf_max is None or ctcf_max < index:
                        ctcf_max = index

            alignments[orientation]['positions'] = [ctcf_min, ctcf_max]

        score_difference = alignments['+']['score'] - alignments['-']['score']

        if abs(score_difference) >= 5:
            if score_difference > 0:
                best = '+'
            else:
                best = '-'
            abs_ctcf_start = alignments[best]['positions'][0] + start
            abs_ctcf_end = alignments[best]['positions'][1] + start

            print(chr, abs_ctcf_start, abs_ctcf_end, '.', record[5], best,
                sep = '\t')
