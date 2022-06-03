#!/bin/bash

: <<'END'
Prepare the bam file from alignment for mutation calling using Mutect2 pipeline
argument(s): bam file (mapping reads to reference), read group name
output(s): sorted and calibrated bam file
END

BAM_FILE=$1
READ_GRP=$2
OUT_DIR=$3
OUT_FILENAME=$4

############ extra parameters ############
basename=${BAM_FILE##*/}
record_name=${basename%.*}

path_to_database="/data/wuchh/somatic-mutation-calls/databases/"
knownSites_Ref1="${path_to_database}1000G_phase1.snps.high_confidence.hg38.vcf"
knownSites_Ref2="${path_to_database}Homo_sapiens_assembly38.dbsnp138.vcf"
intervals="${path_to_database}wgs_calling_regions.hg38.interval_list"
reference="/data/wuchh/somatic-mutation-calls/genomes/hg38/hg38.fa"
##########################################

printf "Processing bam file for: ${record_name}\\n"

## Sort BAM
samtools sort -@ $SLURM_CPUS_PER_TASK -o ${BAM_FILE} ${BAM_FILE}

## AddOrReplaceReadGroups
gatk AddOrReplaceReadGroups -I ${BAM_FILE} -O ${OUT_DIR}${record_name}.rg.bam -RGID ${READ_GRP} -RGLB ${READ_GRP} -RGPL illumina -RGPU ${READ_GRP} -RGSM ${READ_GRP}

## Re-create indicies
samtools index ${OUT_DIR}${record_name}.rg.bam

## MarkDuplicatesSpark
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms60G -Xmx60G -XX:ParallelGCThreads=4" MarkDuplicatesSpark \
-I ${OUT_DIR}${record_name}.rg.bam \
-L ${intervals} \
-O ${OUT_DIR}${record_name}.markdupspark.bam \
-M ${OUT_DIR}${record_name}.metrics.txt
## BaseRecalibrator
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms60G -Xmx60G -XX:ParallelGCThreads=4" BaseRecalibrator \
-I ${OUT_DIR}${record_name}.markdupspark.bam \
-R ${reference} \
--known-sites ${knownSites_Ref1} \
--known-sites ${knownSites_Ref2} \
-L ${intervals} \
-O ${OUT_DIR}${record_name}.recal.table

## ApplyBQSR
gatk --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms60G -Xmx60G -XX:ParallelGCThreads=4" ApplyBQSR \
-R ${reference} \
-I ${OUT_DIR}${record_name}.markdupspark.bam \
--bqsr-recal-file ${OUT_DIR}${record_name}.recal.table \
-L ${intervals} \
-O ${OUT_DIR}${OUT_FILENAME}.recal.bam


