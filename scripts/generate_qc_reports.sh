#!/bin/bash

# USAGE: bash 01_map_pairedendFastQ_bowtie2.sh input/fastq/hgStomachF35 output/hgStomachF35 mm10

# Set user defined environment variables
jobdir="/home/rad/users/gaurav/projects/seqAnalysis/atacseq"
fastqdir=$1
outDir=$2

# Get relevant dirs
bamDir="${outDir}/bams/trimmed"
origFastqcDir="${outDir}/qc/fastqc/00_original"
trimFastqcDir="${outDir}/qc/fastqc/01_trimmed"
bamStatsDir="${outDir}/qc/bamStats"
multiqcDir="${outDir}/qc/multiqc"
mappingLogsDir="${outDir}/logs/trimmingLogs"
trimmedFastDir="${outDir}/tempLocal/trimmed_fastq"

# Create required dirs
echo    "${bamDir} ${mappingLogsDir} ${trimmedFastDir} ${origFastqcDir} ${trimFastqcDir} ${bamStatsDir} ${multiqcDir}"
mkdir -p ${bamDir} ${mappingLogsDir} ${trimmedFastDir} ${origFastqcDir} ${trimFastqcDir} ${bamStatsDir} ${multiqcDir}

# Generate initial quality reports
# Get fastqc reports on original and trimmed fastq files
ls ${fastqdir}/*.fastq.gz | parallel --progress --eta -j 16 "fastqc -o ${origFastqcDir} {}"

# Get multiqc for on original fastqc reports
multiqc -n 01_original_fastq -o ${multiqcDir} ${origFastqcDir}

# Get mapping statistics
for f in $(ls ${bamDir}/*_rmdup.bam); do echo ${f}; samtools stats -@ 64 ${f} > ${bamStatsDir}/$(basename ${f} .bam)_samtools_stats.txt; done;

# Get multiqc for on mapping quality reports
multiqc -n 02_mapping_qc -o ${multiqcDir} ${bamStatsDir}
