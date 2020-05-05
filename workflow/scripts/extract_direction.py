#!/usr/bin/env python3

import os
import re
import sys
import subprocess
import pyCommonTools as pct
from tempfile import NamedTemporaryFile

genome = sys.argv[1]
file = sys.argv[2]


def makeCTCFFasta():
    """ Write 2 temporary FASTA files with forward and reverse-complement
        CTCF consensus sequence.
    """
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5291250/
    ctcf = {'+' : None, '-' : None}
    for ori in ctcf:
        seq = 'CCACNAGGTGGCAG' if ori == '+' else 'CTGCCACCTNGTGG'
        temp = NamedTemporaryFile(delete=False)
        temp.write(f'>ctcf\n{seq}'.encode('utf-8'))
        temp.close()
        ctcf[ori] = temp.name
    return ctcf


def writeFasta(sequence):
    """ Write FASTA sequence string to temporary file. """

    temp = NamedTemporaryFile(delete=False)
    temp.write(sequence)
    temp.close()
    return temp.name


with open(file) as f:

    CTCF = makeCTCFFasta()

    for record in f:
        if record.startswith('#'):
            continue

        record = record.strip().split()
        chr = record[0].strip('chr')
        start = int(record[1])
        end = int(record[2])
        coordinates = f'{chr}:{start}-{end}'

        faidx_command = ['samtools', 'faidx', genome, coordinates]
        faidx = subprocess.Popen(faidx_command,
            stdout = subprocess.PIPE, stderr = subprocess.DEVNULL)

        sequence_path = writeFasta(faidx.communicate()[0])

        if faidx.returncode != 0:
            continue

        alignments = {'+' : {}, '-' : {}}

        for orientation in CTCF:

            command = ['needle', '-gapopen', '10', '-gapextend', '0.5',
                       '-asequence', sequence_path,
                       '-bsequence', CTCF[orientation],
                       '-outfile',  '/dev/stdout']

            needle = subprocess.Popen(
                command, stdout = subprocess.PIPE, stderr = None)

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

            print(chr, abs_ctcf_start, abs_ctcf_end, '.', record[4], best,
                sep = '\t')

    for CTCF_fasta in CTCF.values():
        os.remove(CTCF_fasta)
