#!/bin/bash

# USAGE: bash scripts/parse_nfatac_consensus_peaks_annotation.sh mouse anja nfTALLMm /media/rad/HDD1/atacseq /media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedLibrary/macs/broadPeak/consensus/deseq2/consensus_peaks.mLb.clN.featureCounts.txt /media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.boolean.txt /home/rad/users/gaurav/projects/seqAnalysis/atacseq

species=${1:-"mouse"}
user=${2:-"anja"}
projName=${3:-"rep_tALLcellLineMm"}
outdir=${4:-"/media/rad/HDD1/atacseq"}
origConsFile=${5:-""} # /media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedLibrary/macs/broadPeak/consensus/deseq2/consensus_peaks.mLb.clN.featureCounts.txt
origAnnFile=${6:-""}  # /media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.boolean.txt
jobdir=${7:-"/home/rad/users/gaurav/projects/seqAnalysis/atacseq"}

# Get relevant directories
projDir="${outdir}/${user}/${projName}"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}

# Get parsed consensus peaks bed (Atacseq peaks are broad peaks)
echo "- Getting parsed consensus peaks bed (Atacseq peaks are broad peaks)"
origConsBed="$(dirname ${origAnnFile})/$(basename ${origAnnFile} .boolean.txt).bed"
consensusPeaksBed="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.bed"
awk 'BEGIN { OFS="\t" }{$3=$3 OFS $4"_"$1"_"$2"_"$3}1' ${origConsBed} | cut -f1-4,6,7 | sort -k1,1V -k2,2g > ${consensusPeaksBed}

# Get the tab seaprated raw counts matrix for consensus peaks on all samples
echo "- Getting the tab seaprated raw counts matrix for consensus peaks on all samples)"
rawCountsTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_rawCounts.txt"
sed 1d ${origConsFile} | awk 'BEGIN { OFS="\t" }{$4=$4 OFS $1"_"$2"_"$3"_"$4}1' | cut -f2-5,8- > ${rawCountsTxtFile}

# Sort the input file with the header line intact
(head -n 1 ${rawCountsTxtFile} && tail -n +2 ${rawCountsTxtFile} | sort -k1,1V -k2,2g -k3,3g) > ${rawCountsTxtFile}.tmp && mv ${rawCountsTxtFile}.tmp ${rawCountsTxtFile}

# Remove additional characters from the header
sed -i "s/.mLb.clN.bam//g" ${rawCountsTxtFile}
sed -i "s/.mLb.clN.sorted.bam//g" ${rawCountsTxtFile}

# Rename columns. \b is word boundary to search for full word only
sed -i "s/Chr\b/PeakChrom/" ${rawCountsTxtFile}
sed -i "s/Start\b/PeakStart/" ${rawCountsTxtFile}
sed -i "s/End\b/PeakEnd/" ${rawCountsTxtFile}
sed -i "s/Geneid_Chr_Start_End\b/PeakID/" ${rawCountsTxtFile}
sed -i 's/^/chr/' ${rawCountsTxtFile}
sed -i 's/^/chr/' ${consensusPeaksBed}
sed -i "s/chrPeakChrom\b/PeakChrom/" ${rawCountsTxtFile}
sed -i "s/chrPeakChrom\b/PeakChrom/" ${consensusPeaksBed}

# Get parsed boolean annotation file
peaksAnnTabFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.tab"
Rscript -e 'library(data.table); inputfile  <- commandArgs(TRUE)[1]; outputfile <- commandArgs(TRUE)[2]; peaksDT    <- fread(inputfile, header=TRUE, sep="\t"); dropCols   <- grep("(fc|qval|pval)$", colnames(peaksDT)); peaksDT[,(dropCols) :=NULL]; to.replace <- names(which(sapply(peaksDT, is.logical))); for (var in to.replace) peaksDT[, (var):= as.numeric(get(var))]; 
peaksDT[,PeakID:=paste(interval_id,chr,start,end, sep="_")]; peaksDT[,"interval_id":=NULL]; 
setnames(peaksDT, old = c("chr","start", "end", "num_peaks", "num_samples"), new = c("PeakChrom", "PeakStart", "PeakEnd","NumPeaks","NumSamples")); newColOrder <-  c(colnames(peaksDT)[0:3],"PeakID", colnames(peaksDT)[5:length(colnames(peaksDT))-1]); setcolorder(peaksDT, newColOrder);names(peaksDT) = gsub(pattern = ".mLb.clN.", replacement = "_", x = names(peaksDT));peaksDT[,PeakChrom := paste0("chr",PeakChrom)];fwrite(peaksDT, outputfile, sep = "\t");' ${origAnnFile} ${peaksAnnTabFile}

# Add genomic annotation using chipseeker
peaksAnnTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.txt"

# Create a log file
analysisLogDir="${analysisDir}/logs/analysisLogs"; mkdir -p ${analysisLogDir}
mkdir -p ${analysisLogDir}
annotationLogFile="${analysisLogDir}/$(basename ${consensusPeaksBed} .bed)_annotation.log"

# Annotate consensus peaks
echo "species  = ${species}"  2>&1 | tee    ${annotationLogFile}
echo "user     = ${user}"     2>&1 | tee -a ${annotationLogFile}
echo "projName = ${projName}" 2>&1 | tee -a ${annotationLogFile}
echo "outdir   = ${outdir}"   2>&1 | tee -a ${annotationLogFile}
echo "jobdir   = ${jobdir}"   2>&1 | tee -a ${annotationLogFile}

echo "-------------------------------------------------------------------"
echo "Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksAnnTabFile} -of=${peaksAnnTxtFile} -sp=${species} 2>&1 | tee -a ${annotationLogFile}"
echo "-------------------------------------------------------------------"
Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksAnnTabFile} -of=${peaksAnnTxtFile} -sp=${species} 2>&1 | tee -a ${annotationLogFile}

# Remove the peaksAnnTabFile file as the data is there in the peaksAnnTxtFile
rm ${peaksAnnTabFile}


