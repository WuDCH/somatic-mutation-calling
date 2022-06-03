#!/bin/bash

: <<'END'
Align paired-end reads to the reference genome, hg38 (can be any reference)
argument(s): 2 fastq files, from trimming/quality control
output(s): bam file containing the information, reads mapped to the reference
END

FASTQ1=$1
FASTQ2=$2
OUT_DIR=$3
OUT_FILENAME=$4

#bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /data/wuchh/somatic-mutation-calls/genomes/hg38/hg38.fa ${FASTQ1} ${FASTQ2} > ${OUT_DIR}/${OUT_FILENAME}
bwa mem '/data/wuchh/somatic-mutation-calls/genomes/hg38/hg38.fa' ${FASTQ1} ${FASTQ2} -t $SLURM_CPUS_PER_TASK > ${OUT_DIR}/${OUT_FILENAME}
