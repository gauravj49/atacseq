#!/bin/bash

# USAGE: bash 01_map_pairedendFastQ_bowtie2.sh input/fastq/hgStomachF35 output/hgStomachF35 mm10

species=${1:-"mm10"}
user=${2:-"anja"}
projName=${3:-"rep_tALLcellLineMm"}
outdir=${4:-"/media/rad/HDD1/atacseq"}
jobdir=${5:-"/home/rad/users/gaurav/projects/seqAnalysis/atacseq"}

# Get relevant directories
projDir="${outdir}/${user}/${projName}"
peaksdir="${projDir}/peaks"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
peakFilesDir="${analysisDir}/broadPeaks"; mkdir -p ${peakFilesDir}

# Link broad peaks into a separate directory
ln -s ${peaksdir}/*/macs2peaks/*.broadPeak ${peakFilesDir}

# # Sort the peak files
# for f in ${peakFilesDir}/*.broadPeak; do sort -k1,1 -k2,2n ${f} -o ${f}; done

# Get the peak filenames in a file
ls ${peakFilesDir}/* > ${analysisDir}/peaksFileList.txt

# Get consensus peaks (Atacseq peaks are broad peaks)
echo "- Getting consensus peaks (Atacseq peaks are broad peaks)"
consensusPeaksBed="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.bed"
# mergecols='2,3,4,5,6,7,8,9' # for broad peaks
# collapsecols=$(python -c "import sys; print(','.join(['collapse']*8))")
cat ${peakFilesDir}/*.broadPeak | sort -k1,1 -k2,2n > ${consensusPeaksBed}.tmp
mergeBed -i ${consensusPeaksBed}.tmp > ${consensusPeaksBed}
rm ${consensusPeaksBed}.tmp

# Get the tab seaprated raw counts of the merged peaks for all the samples
# Using deeptools
rawCountsBname="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_rawCounts"
multiBamSummary BED-file --BED ${consensusPeaksBed} --bamfiles ${projDir}/bams/trimmed/*.bam --smartLabels -out ${rawCountsBname}.npz --outRawCounts ${rawCountsBname}.tab -p 64

# Sort the input file with the header line intact
rawCountsTabFile="${rawCountsBname}.tab"
(head -n 1 ${rawCountsTabFile} && tail -n +2 ${rawCountsTabFile} | sort -k1,1V -k2,2g -k3,3g) > ${rawCountsTabFile}.tmp && mv ${rawCountsTabFile}.tmp ${rawCountsTabFile}
