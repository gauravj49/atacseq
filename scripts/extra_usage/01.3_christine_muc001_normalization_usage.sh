# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Normalize raw matrix with vst
R

suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))

inputFile  <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.matrix"
outputFile <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/5320_53631_53646_6075_merge_master_peaks_vstCounts.matrix"

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

