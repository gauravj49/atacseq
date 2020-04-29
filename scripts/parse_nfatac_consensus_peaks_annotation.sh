#!/bin/bash

# USAGE: bash 01_map_pairedendFastQ_bowtie2.sh input/fastq/hgStomachF35 output/hgStomachF35 mm10

species=${1:-"mouse"}
user=${2:-"anja"}
projName=${3:-"rep_tALLcellLineMm"}
outdir=${4:-"/media/rad/HDD1/atacseq"}
origConsFile=${5:-""} # /media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedReplicate/macs/broadPeak/consensus/deseq2/consensus_peaks.mRp.clN.featureCounts.txt
origAnnFile=${6:-""}  # /media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedReplicate/macs/broadPeak/consensus/consensus_peaks.mRp.clN.boolean.annotatePeaks.txt
jobdir=${7:-"/home/rad/users/gaurav/projects/seqAnalysis/atacseq"}

############################
# Remove later
species="mouse"
user="anja"
projName="nfTALLMm"
outdir="/media/rad/HDD1/atacseq"
origConsFile="/media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedReplicate/macs/broadPeak/consensus/deseq2/consensus_peaks.mRp.clN.featureCounts.txt"
origAnnFile="/media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedReplicate/macs/broadPeak/consensus/consensus_peaks.mRp.clN.boolean.txt"
jobdir="/home/rad/users/gaurav/projects/seqAnalysis/atacseq"
############################

# Get relevant directories
projDir="${outdir}/${user}/${projName}"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}

# Get parsed consensus peaks bed (Atacseq peaks are broad peaks)
echo "- Getting parsed consensus peaks bed (Atacseq peaks are broad peaks)"
origConsBed="$(dirname ${origAnnFile})/$(basename ${origAnnFile} .boolean.annotatePeaks.txt).bed"
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

# Rename columns. \b is word boundary to search for full word only
sed -i "s/Chr\b/PeakChrom/" ${rawCountsTxtFile}
sed -i "s/Start\b/PeakStart/" ${rawCountsTxtFile}
sed -i "s/End\b/PeakEnd/" ${rawCountsTxtFile}
sed -i "s/Geneid_Chr_Start_End\b/PeakID/" ${rawCountsTxtFile}

# Get parsed boolean annotation file
peaksAnnTabFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.tab"
Rscript -e 'library(data.table); inputfile  <- commandArgs(TRUE)[1]; outputfile <- commandArgs(TRUE)[2]; peaksDT    <- fread(inputfile, header=TRUE, sep="\t"); dropCols   <- grep("(fc|qval|pval)$", colnames(peaksDT)); peaksDT[,(dropCols) :=NULL]; (to.replace <- names(which(sapply(peaksDT, is.logical)))); for (var in to.replace) peaksDT[, (var):= as.numeric(get(var))]; 
peaksDT[,PeakID:=paste(interval_id,chr,start,end, sep="_")]; setnames(peaksDT, c("PeakChrom", "PeakStart", "PeakEnd","PeakID","NumPeaks","NumSamples"));fwrite(peaksDT, outputfile, sep = "\t");' ${origAnnFile} ${peaksAnnTabFile}


############# Annotation Matrix File #############
# Create binary peaks flag annotation matrix
peaksBedFile="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks_with_peaksID.bed"
peaksAnnTabFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.tab"
tempPkAnnDir="${peakFilesDir}/tempPkAnn"; mkdir -p ${tempPkAnnDir}

# Get the bed file
cut -f1-4 ${rawCountsTxtFile} | sed '1d' > ${peaksBedFile}
colstofilter="," # Column numbers to be extracted at the end
i=1              #
header="PeakChrom\tPeakStart\tPeakEnd\tPeakID\t"
for p in ${peakFilesDir}/*.broadPeak;        
do
 bname=$(basename ${p} _peaks.broadPeak)
 outfile=${tempPkAnnDir}/${bname}.bed
 intersectBed -a ${peaksBedFile} -b ${p} -c | awk '{print $1,$2,$3,$4,$NF}' > ${outfile}
 colstofilter=$(echo "${colstofilter}$((i*5)),");
 header=$(echo -e "${header}${bname}\t")
 echo "${i}) ${bname}: ${colstofilter}"
 echo ""
 i=$((i+1));
done
rm ${peaksBedFile}

# Remove last "," from the from the colstofilter and last tab from the header
colstofilter=$(echo ${colstofilter}|sed 's/,$//')
# ,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90

# Paste output the file with multiple delimiters
# Columns from one files are space separated and multiple files are separated by space
# Using sed to convert the tabs to spaces and then using cut to get the final columns
paste ${tempPkAnnDir}/*.bed| sed 's/\t/ /g' | cut -d ' ' --output-delimiter=$'\t' -f1-4${colstofilter}> ${peaksAnnTabFile}

# Add the header to the file
sed  -i "1i${header}" ${peaksAnnTabFile}

# Removing the last tab in the header
sed -i 's/\t$//' ${peaksAnnTabFile} 

# Add genomic annotation using chipseeker
peaksAnnTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.txt"

# Create a log file
analysisLogDir="${projDir}/logs/analysisLogs"; mkdir -p ${analysisLogDir}
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


