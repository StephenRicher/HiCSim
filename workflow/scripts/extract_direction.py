#!/usr/bin/env python3

import re
import sys
import subprocess
import pyCommonTools as pct


genome = sys.argv[1]
file = sys.argv[2]
chr_ref = sys.argv[3]
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5291250/
ctcf_seqs = {'+' : '>ctcf\nCCACNAGGTGGCAG',
             '-' : '>ctcf\nCTGCCACCTNGTGG'}

with open(file) as f:
    for record in f:
        if record.startswith('#'):
            continue

        record = record.strip().split()
        chr = record[1].strip('chr')
        start = int(record[2])
        end = int(record[3])
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

            # Record index positions of start and end of CTCF alignment
            matches = list(re.finditer(r'[^-]', ctcf_alignment_line))
            ctcf_start = matches[0].start()
            ctcf_end = matches[-1].start()
            alignments[orientation]['positions'] = [ctcf_start, ctcf_end]

        score_difference = alignments['+']['score'] - alignments['-']['score']

        if abs(score_difference) >= 5:
            best = '+' if score_difference > 0 else '-'
            abs_ctcf_start = alignments[best]['positions'][0] + start
            abs_ctcf_end = alignments[best]['positions'][1] + start

            print(chr, abs_ctcf_start, abs_ctcf_end, '.', record[5], best,
                sep = '\t')
