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
multiBamSummary BED-file --BED /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks.bed --bamfiles /media/rad/HDD1/atacseq/sabrina/nkTimecourse/mapping/bams/trimmed/*.bam --smartLabels -out /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.npz --outRawCounts /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab -p 64


# Add peaknames to the file
ipython
#****************************************************************************************************
pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

input_file  = "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab"
output_file = "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.matrix"
peaksDF = pd.read_csv(input_file, sep="\t")
# Fix column names: From #chr to chr and remove <'>
peaksDF.columns = peaksDF.columns.str.strip().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('#', '').str.replace("\'", '').str.replace("_atacseq_mm_se_rmdup", '')
# rename the first column from #chr to chr
peaksDF.rename(columns={ peaksDF.columns[0]: "chr" }, inplace = True)
# Add peaks names to the dataframe
peaksDF.insert (3, "name", ["atacPeak_{0}".format(x) for x in peaksDF.index.tolist()])
peaksDF['peakID'] = peaksDF['chr'].str.cat(peaksDF['start'].apply(str), sep='_').str.cat(peaksDF['end'].apply(str), sep='_').str.cat(peaksDF['name'], sep='_')
# Get column names 
colNames = peaksDF.columns.tolist()
# Move peakID to the front
colNames.insert(0, colNames.pop(colNames.index('peakID')))
# Reorder columns using df.reindex() function
peaksDF = peaksDF.reindex(columns= colNames)
# Drop additional columns
peaksDF.drop(columns=['chr','start','end','name'], inplace=True)
peaksDF.to_csv(output_file, index=False, header=True, sep="\t", float_format='%.0f')
# Cltr+D+D
#****************************************************************************************************

# Normalization
R
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))

inputFile  <- "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.matrix"
outputFile <- "/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_vstCounts.matrix"

# Import data from featureCounts
countdata <- read.table(inputFile, header=TRUE, row.names=1, check.names = FALSE)
coldata   <- data.frame(row.names=colnames(countdata), samples=colnames(countdata))
dds       <- DESeqDataSetFromMatrix(countData=countdata, design=~1, colData=coldata)

vst         <- varianceStabilizingTransformation(dds, blind=TRUE)  
vsd         <- assay(vst)

# # Remove the minimum to get 0s and 0s again
# vsdMin         <- min(vsd)
# vsd            <- vsd - vsdMin

# Convert rownames as feature column to save it in the csv file
vsd           <- cbind(feature = rownames(vsd), vsd)
rownames(vsd) <- 1:nrow(vsd)

# Save the counts file to output csv file
write.table(vsd   , file = outputFile , row.names = F, sep = '\t', quote = F)

