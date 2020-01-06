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

# Add peaknames to the file
ipython
#****************************************************************************************************
pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

input_file  = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.tab"
output_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_rawCounts.matrix"
peaksDF = pd.read_csv(input_file, sep="\t")
# Fix column names: From #chr to chr and remove <'>
peaksDF.columns = peaksDF.columns.str.strip().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('#', '').str.replace("\'", '').str.replace("_r1_001_rmdup", '')
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

