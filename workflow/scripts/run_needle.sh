#!/usr/bin/env bash

genomic_sequence="${1}"
ctcf_sequence="${2}"

needle -asequence <(echo "${genomic_sequence}") \
       -bsequence <(echo "${ctcf_sequence}") \
       -gapopen 10 -gapextend 0.5 \
       -outfile /dev/stdout
