# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Get analysis dirs
mkdir -p /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis

# Get merged peaks
# Modify the peaks list and output file in the megrePeaks.sh file at the begining
# List of summit.bed is in array in the mergePeaks.sh script
bash scripts/mergePeaks.sh 

# Get the tab seaprated raw counts of the merged peaks for all the samples
# Using deeptools
multiBamSummary BED-file --BED /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks.bed --bamfiles /media/rad/HDD1/atacseq/anja/tALLcellLineMm/bams/trimmed/*.bam --smartLabels -out /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.npz --outRawCounts /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.tab -p 64

# Sort the input file with the header line intact
infile="/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.tab"
(head -n 1 ${infile} && tail -n +2 ${infile} | sort -k1,1V -k2,2g -k3,3g) > ${infile}.tmp && mv ${infile}.tmp ${infile}

# Add peaknames to the file
ipython
#****************************************************************************************************
import dask.dataframe
import datatable as dt

import time
pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

input_file  = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.tab"
outtxt_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.txt"
outmat_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.matrix"

# # Read the input file
# start_time = time.time()
# data = dask.dataframe.read_csv(input_file,sep="\t", header=0)
# print("%s seconds" % (time.time() - start_time))
# # 0.0748898983001709 seconds

# Header starting with # is a problem
# start_time = time.time()
# peaksDT = dt.fread(input_file, sep="\t", header=True, nthreads=16)
# print("%s seconds" % (time.time() - start_time))
# # 0.0748898983001709 seconds

start_time = time.time()
peaksDF = pd.read_csv(input_file, sep="\t")
print("%s seconds" % (time.time() - start_time))
# 6.941580533981323 seconds

# Fix column names: From #chr to chr and remove <'>
peaksDF.columns = peaksDF.columns.str.strip().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('#', '').str.replace("\'", '').str.replace("_atacseq_mm_se_rmdup", '')
# rename the first three column from #chr to PeakChrom, start to PeakStart and end to PeakEnd
peaksDF.rename(columns={ peaksDF.columns[0]: "PeakChrom" }, inplace = True)
peaksDF.rename(columns={ peaksDF.columns[1]: "PeakStart" }, inplace = True)
peaksDF.rename(columns={ peaksDF.columns[2]: "PeakEnd"   }, inplace = True)
# Add peaks names to the dataframe
# peaksDF.insert (3, "name", ["atacPeak_{0}".format(x) for x in peaksDF.index.tolist()])
# peaksDF['peakID'] = peaksDF['PeakChrom'].str.cat(peaksDF['PeakStart'].apply(str), sep='_').str.cat(peaksDF['PeakEnd'].apply(str), sep='_').str.cat(peaksDF['name'], sep='_')
peaksDF.insert (3, "PeakID", peaksDF['PeakChrom'].str.cat(peaksDF['PeakStart'].apply(str), sep='_').str.cat(peaksDF['PeakEnd'].apply(str), sep='_').str.cat(["atacPeak_{0}".format(x) for x in peaksDF.index.tolist()], sep='_'))
# Save the text file
peaksDF.to_csv(outtxt_file, index=False, header=True, sep="\t", float_format='%.0f')

# Drop additional columns
peaksDF.drop(columns=['PeakChrom','PeakStart','PeakEnd'], inplace=True)
peaksDF.to_csv(outmat_file, index=False, header=True, sep="\t", float_format='%.0f')

# Cltr+D+D
#****************************************************************************************************

############# Annotation Matrix File #############
# Create binary peaks flag annotation matrix
peaksTxtFile='/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.txt'
peaksBedFile='/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks.bed'
peaksAnnFile='/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_annotation.tab'
peakFilesDir='/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks'
tempPkAnnDir="${peakFilesDir}/tempPkAnn"; mkdir -p ${tempPkAnnDir}

# #  Copy peaks file
# cp -rv /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/peaks/*/macs2peaks/*.broadPeak ${peakFilesDir}

# Get the bed file
cut -f1-4 ${peaksTxtFile} | sed '1d' > ${peaksBedFile}
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

# Remove last "," from the from the colstofilter and last tab from the header
colstofilter=$(echo ${colstofilter}|sed 's/,$//')
# ,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90

# Paste output the file with multiple delimiters
# Columns from one files are space separated and multiple files are separated by space
# Using sed to convert the tabs to spaces and then using cut to get the final columns
paste ${tempPkAnnDir}/*.bed| sed 's/\t/ /g' | cut -d ' ' --output-delimiter=$'\t' -f1-4${colstofilter}> ${peaksAnnFile}

# Add the header to the file
sed  -i "1i${header}" ${peaksAnnFile}

# Removing the last tab in the header
sed -i 's/\t$//' ${peaksAnnFile} 

# Add genomic annotation using annotater
R
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("annotatr")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.ensGene))
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

# Annotate each sample peaks with ChIPseeker
# Read bed file
cat("\n\t- Reading the bed file\n")
inputfile  <- '/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_annotation.tab'
outputfile <- '/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_annotation.txt'

origPeakDT <- fread(inputfile, header=TRUE, sep="\t")
peaksBed   <- makeGRangesFromDataFrame(as.data.frame(origPeakDT))

origPeakDT <- fread(inputfile, header=TRUE, sep="\t")
peaksBed   <- makeGRangesFromDataFrame(as.data.frame(origPeakDT))
# peaksBed <- readGeneric(origPeakDT,chr=1,start=2,end=3,strand=NULL)

# Annotate regions
cat("\n\t- Annotate regions\n")
peakAnno <- annotatePeak(peaksBed, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# Convert to dataframe 
peakAnnoDT <- as.data.table(peakAnno)

# Rename the header
# Original header: > colnames(as.data.table(peakAnno))
# [1]  "seqnames"      "start"         "end"           "width"        
# [5]  "strand"        "annotation"    "geneChr"       "geneStart"    
# [9]  "geneEnd"       "geneLength"    "geneStrand"    "geneId"       
# [13] "transcriptId"  "distanceToTSS" "ENTREZID"      "SYMBOL"       
# [17] "GENENAME" 
setnames(peakAnnoDT, c("PeakChrom", "PeakStart", "PeakEnd", "PeakWidth", "PeakStrand", "DetailedGenomicAnnotation", "GeneChrom", "GeneStart", "GeneEnd", "GeneLength", "GeneStrand", "GeneID", "TranscriptID", "DistanceToTSS", "EntrezID", "GeneName", "GeneDesc"))

# Copy the DetailedGenomicAnnotation column as GenomicAnnotation column
peakAnnoDT[,GenomicAnnotation:=DetailedGenomicAnnotation]

# Replace the detailed annotation to the abstract annotation
peakAnnoDT[DetailedGenomicAnnotation %like%   'exon 1 '   , GenomicAnnotation:='ExonFirst']
peakAnnoDT[!(DetailedGenomicAnnotation %like% 'exon 1 ')  , GenomicAnnotation:='ExonOther']
peakAnnoDT[DetailedGenomicAnnotation %like%   'intron 1 ' , GenomicAnnotation:='IntronFirst']
peakAnnoDT[!(DetailedGenomicAnnotation %like% 'intron 1 '), GenomicAnnotation:='IntronOther']
peakAnnoDT[DetailedGenomicAnnotation=='Distal Intergenic' , GenomicAnnotation:='IntergenicDistal']
peakAnnoDT[DetailedGenomicAnnotation=="3' UTR"            , GenomicAnnotation:='ThreeUTR']
peakAnnoDT[DetailedGenomicAnnotation=="5' UTR"            , GenomicAnnotation:='FiveUTR' ]
peakAnnoDT[DetailedGenomicAnnotation=='Downstream (1-2kb)', GenomicAnnotation:='DownstreamProximal']
peakAnnoDT[DetailedGenomicAnnotation=='Downstream (<1kb)' , GenomicAnnotation:='DownstreamBasal']
peakAnnoDT[DetailedGenomicAnnotation=='Downstream (2-3kb)', GenomicAnnotation:='DownstreamDistal']
peakAnnoDT[DetailedGenomicAnnotation=='Promoter (1-2kb)'  , GenomicAnnotation:='PromoterProximal']
peakAnnoDT[DetailedGenomicAnnotation=='Promoter (<=1kb)'  , GenomicAnnotation:='PromoterBasal']
peakAnnoDT[DetailedGenomicAnnotation=='Promoter (2-3kb)'  , GenomicAnnotation:='PromoterDistal']

# Reorder the columns
setcolorder(peakAnnoDT, c("PeakChrom", "PeakStart", "PeakEnd", "PeakWidth", "PeakStrand", "GenomicAnnotation", "DetailedGenomicAnnotation", "GeneChrom", "GeneStart", "GeneEnd", "GeneLength", "GeneStrand", "GeneID", "TranscriptID", "DistanceToTSS", "EntrezID", "GeneName", "GeneDesc"))

# peakAnnoDT[,unique(GenomicAnnotation)]
#  [1] "IntronFirst"        "IntronOther"        "PromoterBasal"     
#  [4] "ThreeUTR"           "IntergenicDistal"   "PromoterDistal"    
#  [7] "PromoterProximal"   "DownstreamProximal" "DownstreamBasal"   
# [10] "DownstreamDistal"   "FiveUTR"

# Merge the orginal data table with annotation data table
# There are few entires in the annotation data table that ...
# ... were not present but the output datatable should be of same size as input
# Create Temporary ID for merging
origPeakDT[,mergeID:=paste0(PeakChrom,PeakStart,PeakEnd)] 
peakAnnoDT[,mergeID:=paste0(PeakChrom,PeakStart,PeakEnd)] 

# Merge on common ids and keep all the entires from the first data table
mergedPeaksDT <- merge(peakAnnoDT, origPeakDT,all.y=T)

# Remove the 
mergedPeaksDT[, c('mergeID') :=NULL]

# Move PeakID column to the 4th position
setcolorder(mergedPeaksDT, c("PeakChrom", "PeakStart", "PeakEnd", "PeakID"))

# Save results in the output file
fwrite(mergedPeaksDT, outputfile, sep = "\t")

# Sort alphanumerically with PeakID
system(paste0("sort -k1,1V -k2,2g -k3,3g ", outputfile, " -o ", outputfile))

cat(paste0("\t- ",outputfile,"\n"))

#########################################################################

#****************************************************************************************************

# Normalize raw matrix with vst
R
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
inputFile  <- "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.matrix"
outputFile <- "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_vstCounts.matrix"

# Import data from featureCounts
countdata <- read.table(inputFile, header=TRUE, row.names=1, check.names = FALSE)
coldata   <- data.frame(row.names=colnames(countdata), samples=colnames(countdata))
dds       <- DESeqDataSetFromMatrix(countData=countdata, design=~1, colData=coldata)

vst         <- varianceStabilizingTransformation(dds, blind=TRUE)  
vsd         <- assay(vst)

# Convert rownames as feature column to save it in the csv file
vsd           <- cbind(feature = rownames(vsd), vsd)
rownames(vsd) <- 1:nrow(vsd)

# Save the counts file to output csv file
write.table(vsd   , file = outputFile , row.names = F, sep = '\t', quote = F)

