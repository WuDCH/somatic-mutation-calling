#!/bin/bash

: <<'END'
Retrieve the raw fastq files for paired-end read experiments
argument(s): SRA run accession
output(s): 2 fastq files containing the reads
END

SRR_ACCESSION=$1
OUT_DIR=$2
NGC_FILE=$3

fasterq-dump --ngc ${NGC_FILE} -t /lscratch/$SLURM_JOBID -O ${OUT_DIR} ${SRR_ACCESSION}
