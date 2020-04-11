# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

jobdir="/home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="mm10"
user="anja"
projName="rep_tALLcellLineMm"

# Get relevant directories
outdir="/media/rad/HDD1/atacseq"
projDir="${outdir}/${user}/${projName}"
peaksdir="${projDir}/peaks"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
peakFilesDir="${analysisDir}/broadPeaks"; mkdir -p ${peakFilesDir}

# Get consensus peaks, raw counts and annotation matrix for consensus peaks
bash scripts/get_peaks_counts_annotation_files.sh ${species} ${user} ${projName} ${outdir} ${jobdir}

#****************************************************************************************************
# Normalize raw matrix with vst
R
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("data.table", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("splitstackshape", warn.conflicts=FALSE, quietly=TRUE))

inputFile  <- "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/tALLcellLineMm_mergedreps_all_merge_master_peaks_rawCounts.matrix"
outputFile <- "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/tALLcellLineMm_mergedreps_all_merge_master_peaks_vstCounts.matrix"

# Import data from featureCounts
countdata  <- read.table(inputFile, header=TRUE, row.names=1, check.names = FALSE)

# Create the sampletable from the column names
targetsDT <- data.table(colnames(countdata))
SampleTable <- cSplit(targetsDT, 'V1','_', drop=F)

# Create new column batch column
batches    <- c(1,1,2,2,1,1,2,2,3,3,2,2,2,1)
SampleTable$batches <- batches

# Rename V_1 column to sampleNames
names(SampleTable)[names(SampleTable) == "V1"]   <- "sampleNames"

# Drop useless columns
SampleTable <- subset(SampleTable, select = -c(V1_1, V1_2, V1_3, V1_4, V1_5, V1_6, V1_7, V1_8, V1_9))

# Get the dds object
# coldata   <- data.frame(row.names=colnames(countdata), samples=colnames(countdata))
# dds       <- DESeqDataSetFromMatrix(countData=countdata, design=~1, colData=coldata)

# the design formula contains a numeric variable with integer values,
# specifying a model with increasing fold change for higher values.
# did you mean for this to be a factor? if so, first convert
# this variable to a factor using the factor() function
SampleTable$batches <- as.factor(SampleTable$batches)
dds <- DESeqDataSetFromMatrix(colData  = SampleTable, countData=countdata, design = ~batches)

# VST 
vst         <- varianceStabilizingTransformation(dds, blind=TRUE)  
vsd         <- assay(vst)

# Convert rownames as feature column to save it in the csv file
vsd           <- cbind(feature = rownames(vsd), vsd)
rownames(vsd) <- 1:nrow(vsd)

# Save the counts file to output csv file
write.table(vsd   , file = outputFile , row.names = F, sep = '\t', quote = F)

