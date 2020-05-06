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

inputFile  <- "/media/rad/HDD1/atacseq/anja/rep_tALLcellLineMm/analysis/rep_tALLcellLineMm_analysis_consensus_peaks_rawCounts.txt"
outputFile <- "/media/rad/HDD1/atacseq/anja/rep_tALLcellLineMm/analysis/rep_tALLcellLineMm_analysis_consensus_peaks_vstCounts.txt"

# Import data from featureCounts
countdata  <- fread(inputFile, header=TRUE)

# Subset genomic annoations from the counts DT
annotDT <- countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID')]

# Remove PeakChrom PeakStart  PeakEnd columns
countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID'):=NULL]  

# Create the sampletable from the column names
targetsDT <- data.table(colnames(countdata))
SampleTable <- cSplit(targetsDT, 'V1','_', drop=F)

# Create new column batch column
batches    <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1,1,2,2,1,1,1,1,1,1)
SampleTable$batches <- batches

# Rename V1_1 column to sampleNames
names(SampleTable)[names(SampleTable) == "V1_1"]   <- "sampleNames"

# Drop useless columns
SampleTable <- subset(SampleTable, select = -c(V1_2, V1_3, V1_4, V1_5, V1_6))

# Get the dds object
# the design formula contains a numeric variable with integer values,
# specifying a model with increasing fold change for higher values.
# did you mean for this to be a factor? if so, first convert
# this variable to a factor using the factor() function
SampleTable$batches <- as.factor(SampleTable$batches)
dds <- DESeqDataSetFromMatrix(colData  = SampleTable, countData=countdata, design = ~batches)

# VST 
vst      <- varianceStabilizingTransformation(dds, blind=TRUE)  
vsd      <- assay(vst)

# Add the genomic annotation in base R
vstBindDT <- cbind(annotDT, vsd)

# Round all normalized counts to 4 decimal places
vstOutDT  <- data.table(vstBindDT %>% mutate_at(vars(-PeakChrom,-PeakStart,-PeakEnd,-PeakID), funs(round(., 4))))

# Save the counts file to output csv file
fwrite(vstOutDT   , file = outputFile , sep = '\t', quote = F)

