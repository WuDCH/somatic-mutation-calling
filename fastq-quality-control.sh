#!/bin/bash

: <<'END'
Perform quality control, remove or trim read/read segments/adaptors
given default thresholds. On paired-end reads
argument(s): 2 fastq files, *_1.fastq and *_2.fastq
output(s): 2 fastq files, containing the trimmed reads
END

FASTQ1=$1
FASTQ2=$2
OUT_DIR=$3

trim_galore --paired --fastqc --fastqc_args "--outdir ${OUT_DIR}" --cores 4 ${FASTQ1} ${FASTQ2} -o ${OUT_DIR} --no_report_file
rm ${FASTQ1} ${FASTQ2} # delete raw fastq files to free storage
