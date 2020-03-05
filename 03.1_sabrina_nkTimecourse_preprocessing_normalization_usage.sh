# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# One sample was named wrong 
mv GS400_71h_NKT_GGACTCCT_atacseq_mm_se.fastq.gz GS400_72h_NKT_GGACTCCT_atacseq_mm_se.fastq.gz

# Batch3 processing
mv B3_L3_TAAGGCGA_R1.fastq.gz SB01_DP_TAAGGCGA_atacseq_mm_se.fastq.gz
mv B3_L3_CGAGGCTG_R1.fastq.gz SB03_24h_NKT_CGAGGCTG_atacseq_mm_se.fastq.gz
mv B3_L4_CGAGGCTG_R1.fastq.gz SB04_24h_NKT_CGAGGCTG_atacseq_mm_se.fastq.gz
mv B3_L3_AAGAGGCA_R1.fastq.gz SB06_36h_NKT_AAGAGGCA_atacseq_mm_se.fastq.gz
mv B3_L4_GTAGAGGA_R1.fastq.gz SB11_48h_NKT_DNAse_GTAGAGGA_atacseq_mm_se.fastq.gz
mv B3_L3_GTCGTGAT_R1.fastq.gz SB08_36h_NKT_DNAse_GTCGTGAT_atacseq_mm_se.fastq.gz
mv B3_L4_GTCGTGAT_R1.fastq.gz SB17_72h_NKT_GTCGTGAT_atacseq_mm_se.fastq.gz
mv B3_L4_ACCACTGT_R1.fastq.gz SB19_72h_NKT_DNAse_ACCACTGT_atacseq_mm_se.fastq.gz
mv B3_L3_TGGATCTG_R1.fastq.gz SB05_24h_NKT_DNAse_TGGATCTG_atacseq_mm_se.fastq.gz
mv B3_L3_CCGTTTGT_R1.fastq.gz SB07_36h_NKT_DNAse_CCGTTTGT_atacseq_mm_se.fastq.gz
mv B3_L3_CGTACTAG_R1.fastq.gz SB09_48h_NKT_CGTACTAG_atacseq_mm_se.fastq.gz
mv B3_L4_CGTACTAG_R1.fastq.gz SB10_48h_NKT_CGTACTAG_atacseq_mm_se.fastq.gz
mv B3_L5_CGTACTAG_R1.fastq.gz SB12_60h_NKT_CGTACTAG_atacseq_mm_se.fastq.gz
mv B3_L3_AGGCAGAA_R1.fastq.gz SB13_60h_NKT_AGGCAGAA_atacseq_mm_se.fastq.gz
mv B3_L4_AGGCAGAA_R1.fastq.gz SB14_60h_NKT_DNAse_AGGCAGAA_atacseq_mm_se.fastq.gz
mv B3_L4_TCCTGAGC_R1.fastq.gz SB15_72h_NKT_TCCTGAGC_atacseq_mm_se.fastq.gz
mv B3_L4_GGACTCCT_R1.fastq.gz SB21_84h_NKT_GGACTCCT_atacseq_mm_se.fastq.gz
mv B3_L4_TAGGCATG_R1.fastq.gz SB22_96h_NKT_TAGGCATG_atacseq_mm_se.fastq.gz
mv B3_L3_CTCTCTAC_R1.fastq.gz SB16_72h_NKT_CTCTCTAC_atacseq_mm_se.fastq.gz
mv B3_L4_CTCTCTAC_R1.fastq.gz SB18_72h_NKT_CTCTCTAC_atacseq_mm_se.fastq.gz
mv B3_L5_CTCTCTAC_R1.fastq.gz SB20_72h_NKT_DNAse_CTCTCTAC_atacseq_mm_se.fastq.gz
mv B3_L3_CAGAGAGG_R1.fastq.gz SB02_DP_CAGAGAGG_atacseq_mm_se.fastq.gz

# 1) Generate initial quality reports
cd /media/rad/HDD1/atacseq/sabrina/nkTimecourse
mkdir -p qc/{fastqc,bamStats}

# Get fastqc reports on original and trimmed fastq files
ls fastq/*.fastq.gz | parallel --progress --eta -j 16 'fastqc -o qc/fastqc/00_original {}'
ls mapping/tempLocal/trimmed_fastq/*.fq.gz | parallel --progress --eta -j 16 'fastqc -o qc/fastqc/01_trimmed {}'

# Get multiqc for on original fastqc reports
multiqc -n 01_original_fastq -o qc/multiqc qc/fastqc/00_original

# Get mapping statistics
for f in $(ls mapping/bams/trimmed/*.bam); do echo ${f}; samtools stats -@ 64 ${f} > qc/bamStats/$(basename ${f} .bam)_samtools_stats.txt; done;

# Get multiqc for on mapping quality reports
multiqc -n 02_mapping_qc -o qc/multiqc qc/bamStats

# Get merged peaks
# Get analysis dirs
mkdir -p /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis

# modify the peaks list and output file in the megrePeaks.sh file at the begining
# List of summit.bed is in array in the mergePeaks.sh script
bash scripts/mergePeaks.sh 

# Get the tab seaprated raw counts of the merged peaks for all the samples
# Using deeptools
multiBamSummary BED-file --BED /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks.bed --bamfiles /media/rad/HDD1/atacseq/sabrina/nkTimecourse/bams/trimmed/*.bam --smartLabels -out /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.npz --outRawCounts /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab -p 64

# Sort the input file 
# Make sure to check the header during sorting
sort -k1,1V -k2,2g -k3,3g "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab" -o "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab"

# Add peaknames to the file
ipython
#****************************************************************************************************
pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

input_file  = "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab"
outtxt_file = "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.txt"
outtxt_file = "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.txt"

peaksDF = pd.read_csv(input_file, sep="\t")
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
peaksTabFile="/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.txt"
peaksBedFile="/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks.bed"
peaksAnnFile='/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_annotation.tab'
peakFilesDir='/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/peakFiles'
tempPkAnnDir="${peakFilesDir}/tempPkAnn"; mkdir -p ${tempPkAnnDir}

# Get the bed file
cut -f1-4 ${peaksTabFile} | sed '1d' > ${peaksBedFile}
colstofilter="," # Column numbers to be extracted at the end
i=1              #
header="PeakChrom\tPeakStart\tPeakEnd\tPeakID\t"        
for p in ${peakFilesDir}/*.broadPeak;
do
  bname=$(basename ${p} _atacseq_mm_se_rmdup_peaks.broadPeak)
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
paste /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/peakFiles/tempPkAnn/*.bed| sed 's/\t/ /g' | cut -d ' ' --output-delimiter=$'\t' -f1-4${colstofilter}> ${peaksAnnFile}

# Add the header to the file
sed  -i "1i${header}" ${peaksAnnFile}

# Removing the last tab in the header
sed -i 's/\t$//' ${peaksAnnFile} 

# # Get the bed file without the peaks column
# peaksAnnBed='/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_annotation.bed'
# cut -f1-4 ${peaksAnnFile} > ${peaksAnnBed}
# # Remove header
# sed -i '1d' ${peaksAnnBed}

# Add genomic annotation using chipseeker
R
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("annotatr")
# suppressPackageStartupMessages(library(annotatr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.ensGene))
# suppressPackageStartupMessages(library(ReactomePA))
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene 

# Annotate each sample peaks with ChIPseeker
# Read bed file
cat("\n\t- Reading the bed file\n")
inputfile  <- '/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_annotation.tab'
outputfile <- '/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_annotation.txt'

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
peakAnnoDT[DetailedGenomicAnnotation=='Downstream (1-2kb)', GenomicAnnotation:='DownstreamBasal']
peakAnnoDT[DetailedGenomicAnnotation=='Downstream (<1kb)' , GenomicAnnotation:='DownstreamProximal']
peakAnnoDT[DetailedGenomicAnnotation=='Downstream (2-3kb)', GenomicAnnotation:='DownstreamDistal']
peakAnnoDT[DetailedGenomicAnnotation=='Promoter (1-2kb)'  , GenomicAnnotation:='PromoterBasal']
peakAnnoDT[DetailedGenomicAnnotation=='Promoter (<=1kb)'  , GenomicAnnotation:='PromoterProximal']
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

# Normalization
# # Theory:
# "The rlog transformation is calculated by fitting for each gene a GLM with a baseline expression (i.e., intercept only) and, computing for each sample, shrunken LFCs with respect to the baseline, using the same empirical Bayes procedure as before (Materials and methods)." 
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8. So instead of using the design, the rlog creates a matrix with an intercept and a coefficient for each sample. See Methods of DESeq2 paper for more details.

# Whereas the VST is conceptually different. It looks at the trend between variance and mean in the data, and then tries to find a strictly monotonous transformation of the data so that this trend is removed. In practice, the transformation will approach the logarithm function for high values and the square root function for small values (incl. 0), and smoothly interpolate inbetween.

# Source: https://support.bioconductor.org/p/104615/

R
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("data.table", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("splitstackshape", warn.conflicts=FALSE, quietly=TRUE))

inputFile  <- "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.matrix"

# Import data from featureCounts
countdata <- read.table(inputFile, header=TRUE, row.names=1, check.names = FALSE)

# Replace sample name from 71 to 72 as it belongs to 72h timepoint
colnames(countdata) <- gsub('71','72', colnames(countdata))

# Create the sampletable from the column names
targetsDT <- data.table(colnames(countdata))
SampleTable <- cSplit(targetsDT, 'V1','_', drop=F)

# Create new column batch column
# SampleTable$<-gsub('[0-9]+', '', SampleTable$V1_1)
batches <- c(rep('1',8),rep('2',10),rep('3',15))
SampleTable$batches <- batches

# Rename V_2 column to timepoint
names(SampleTable)[names(SampleTable) == "V1_2"] <- "timepoints"
names(SampleTable)[names(SampleTable) == "V1"]   <- "sampleNames"

# Drop useless columns
SampleTable <- subset(SampleTable, select = -c(V1_1, V1_3, V1_4))

# Get the dds object
dds <- DESeqDataSetFromMatrix(colData  = SampleTable, countData=countdata, design = ~timepoints)

# Relative Log Transformation
# Transform data to log space and visualize samples
rld <- rlogTransformation(dds, blind = TRUE)

# Convert rownames as feature column to save it in the csv file
rld           <- cbind(feature = rownames(rld), rld)
rownames(rld) <- 1:nrow(rld)
# Save the counts file to output csv file
outputFile <- "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rlogCounts.matrix"
write.table(rld   , file = outputFile , row.names = F, sep = '\t', quote = F)

# VST 
vst         <- varianceStabilizingTransformation(dds, blind=TRUE)  
vsd         <- assay(vst)
# Convert rownames as feature column to save it in the csv file
vsd           <- cbind(feature = rownames(vsd), vsd)
rownames(vsd) <- 1:nrow(vsd)
# Save the counts file to output csv file
outputFile <- "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_vstCounts.matrix"
write.table(vsd   , file = outputFile , row.names = F, sep = '\t', quote = F)


# Source: https://f1000research.com/articles/5-1408
# Law, C. W., Alhamdoosh, M., Su, S., Dong, X., Tian, L., Smyth, G. K., & Ritchie, M. E. (2016). RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR. F1000Research, 5.

# The density of log-CPM values for raw pre-filtered data (A) and post-filtered data (B) are shown for each sample. Dotted vertical lines mark the log-CPM threshold (equivalent to a CPM value of about 0.2) used in the filtering step.
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}