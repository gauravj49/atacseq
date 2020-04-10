#!/bin/bash

# USAGE: bash 01_map_pairedendFastQ_bowtie2.sh input/fastq/hgStomachF35 output/hgStomachF35 mm10

species=${1:-"mouse"}
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

# # Get peaks summit for broad peaks
# for broadpeak in ${peakFilesDir}/*.broadPeak; 
# do 
#  bname=$(basename ${broadpeak} _peaks.broadPeak); 
#  echo ${bname}; 
#  bamfile=${projDir}/bams/trimmed/${bname}.bam; 
#  summitfile=${peaksdir}/${bname}/macs2peaks/${bname}_peaks.summit
#  macs2 refinepeak -b ${broadpeak} -i ${bamfile} -o ${summitfile}
#  echo ""
# done

# # Link broad peaks into a separate directory
# ln -s ${peaksdir}/*/macs2peaks/*.summit ${peakFilesDir}

# Get the peak filenames in a file
ls ${peakFilesDir}/* > ${analysisDir}/peaksFileList.txt

# Get consensus peaks (Atacseq peaks are broad peaks)
echo "- Getting consensus peaks (Atacseq peaks are broad peaks)"
consensusPeaksBed="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.bed"
cat ${peakFilesDir}/*.broadPeak | sort -k1,1V -k2,2g > ${consensusPeaksBed}.tmp
mergeBed -i ${consensusPeaksBed}.tmp > ${consensusPeaksBed}
rm ${consensusPeaksBed}.tmp

# Get the tab seaprated raw counts of the merged peaks for all the samples
# Using deeptools
rawCountsBname="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_rawCounts"
multiBamSummary BED-file --BED ${consensusPeaksBed} --bamfiles ${projDir}/bams/trimmed/*.bam --smartLabels -out ${rawCountsBname}.npz --outRawCounts ${rawCountsBname}.txt -p 64

# Sort the input file with the header line intact
rawCountsTxtFile="${rawCountsBname}.txt"
(head -n 1 ${rawCountsTxtFile} && tail -n +2 ${rawCountsTxtFile} | sort -k1,1V -k2,2g -k3,3g) > ${rawCountsTxtFile}.tmp && mv ${rawCountsTxtFile}.tmp ${rawCountsTxtFile}

# Remove single quotes from the header
sed -i "s/[']//g" ${rawCountsTxtFile}

# Rename columns. \b is word boundary to search for full word only
sed -i "s/#chr\b/PeakChrom/" ${rawCountsTxtFile}
sed -i "s/start\b/PeakStart/" ${rawCountsTxtFile}
sed -i "s/end\b/PeakEnd/" ${rawCountsTxtFile}

# Change float to integer in the python one liner
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; outtxt_file=sys.argv[2]; peaksDT = dt.fread(input_file, sep='\t', header=True, nthreads=16);peaksDF = peaksDT.to_pandas(); peaksDF.to_csv(outtxt_file, header=True, index=False, sep='\t', float_format='%.0f')" ${rawCountsTxtFile} ${rawCountsTxtFile}.tmp

# Add PeakID column
# Source for concatenating: https://unix.stackexchange.com/questions/304499/adding-column-to-a-table-by-concatenating-values-from-other-columns
# source for line numbers: https://stackoverflow.com/questions/20752043/print-line-numbers-starting-at-zero-using-awk
awk 'BEGIN { OFS="\t" }{$3=$3 OFS "atacPeak_"i++"_"$1"_"$2"_"$3}1' ${rawCountsTxtFile}.tmp | awk 'BEGIN { OFS="\t" }NR==1{$4="PeakID"}1' > ${rawCountsTxtFile} && rm ${rawCountsTxtFile}.tmp

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



