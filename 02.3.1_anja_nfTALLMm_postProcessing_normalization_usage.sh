# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Parameters for the script
species="mouse"
user="anja"
projName="nfTALLMm"
outdir="/media/rad/HDD1/atacseq"
origConsFile="/media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedLibrary/macs/broadPeak/consensus/deseq2/consensus_peaks.mLb.clN.featureCounts.txt"
origAnnFile="/media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.boolean.txt"
jobdir="/home/rad/users/gaurav/projects/seqAnalysis/atacseq"

# Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
consensusPeaksBed="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.bed"
rawCountsTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_rawCounts.txt"
peaksAnnTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.txt"

# 1) Parse concensus raw matrix and boolean matrix to get annoation files
echo "bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}"
bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}

# Output files are:
echo "- Consensus bed file: ${consensusPeaksBed}"
echo "- Raw peaks count   : ${rawCountsTxtFile}"
echo "- Peaks annotation  : ${peaksAnnTxtFile}"

# 2) Normalize raw matrix with vst
R
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("data.table", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("splitstackshape", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("dplyr", warn.conflicts=FALSE, quietly=TRUE))

inputFile  <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_analysis_consensus_peaks_rawCounts.txt"
outputFile <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_analysis_consensus_peaks_vstCounts.txt"

# Import data from featureCounts
countdata  <- fread(inputFile, header=TRUE)

# Subset genomic annoations from the counts DT
annotDT <- countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID')]

# Remove PeakChrom PeakStart  PeakEnd columns
countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID'):=NULL]  

# Create the sampletable from the column names
targetsDT <- data.table(colnames(countdata))
# SampleTable <- cSplit(targetsDT, 'V1','_', drop=F)
SampleTable <- targetsDT

# Create new column batch column
batches    <- c(2,1,1,2,2,2,1,2,1,2,3,1,1,1,2,1,1,2,3,1,1,2,1,2,1,1)
SampleTable$batches <- batches

# Rename V1 column to sampleNames
names(SampleTable)[names(SampleTable) == "V1"]   <- "sampleNames"

# # Drop useless columns
# SampleTable <- subset(SampleTable, select = -c(V1_2, V1_3, V1_4, V1_5, V1_6))

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

###############################################################
#                       DOCUMENTATION
###############################################################

# 1) Parse concensus raw matrix and boolean matrix to get annoation files
# 1.1) Raw matrix header includes

# # Program:featureCounts v1.6.0; Command:"featureCounts" "-F" "SAF" "-O" "--fracOverlap" "0.2" "-T" "6" "-p" "--donotsort" "-a" "consensus_peaks.mLb.clN.saf" "-o" "consensus_peaks.mLb.clN.featureCounts.txt" "DN3Tcell_R2.mLb.clN.bam" "CD8SinglePositive_R1.mLb.clN.bam" "CD4SinglePositive_R1.mLb.clN.bam" "DN2aTcell_R2.mLb.clN.bam" "DN4Tcell_R2.mLb.clN.bam" "DN3Tcell_R1.mLb.clN.bam" "HSCBoneMarrow_R2.mLb.clN.bam" "DN2bTcell_R2.mLb.clN.bam" "BcellBoneMarrow_R2.mLb.clN.bam" "DN4Tcell_R1.mLb.clN.bam" "DPTcell_R1.mLb.clN.bam" "HSCBoneMarrow_R1.mLb.clN.bam" "MPPBoneMarrow_R1.mLb.clN.bam" "NKBoneMarrow_R2.mLb.clN.bam" "ETPTcell_R1.mLb.clN.bam" "BcellBoneMarrow_R1.mLb.clN.bam" "CD4SinglePositive_R2.mLb.clN.bam" "DN2bTcell_R1.mLb.clN.bam" "DPTcell_R2.mLb.clN.bam" "CLPBoneMarrow_R2.mLb.clN.bam" "CLPBoneMarrow_R1.mLb.clN.bam" "DN2aTcell_R1.mLb.clN.bam" "MPPBoneMarrow_R2.mLb.clN.bam" "ETPTcell_R2.mLb.clN.bam" "NKBoneMarrow_R1.mLb.clN.bam" "CD8SinglePositive_R2.mLb.clN.bam"
#  ------------ ----- --------- --------- -------- -------- ------------------------- ---------------------------------- --- 
#     Geneid     Chr    Start      End     Strand   Length   DN3Tcell_R2.mLb.clN.bam   CD4SinglePositive_R1.mLb.clN.bam   …  
#  ------------ ----- --------- --------- -------- -------- ------------------------- ---------------------------------- --- 
 
#   Interval_1     1   3008719   3009020   +           302                         1                                  9   …  
#   Interval_2     1   3042821   3043241   +           421                         7                                  0   …  
#   Interval_3     1   3174534   3174707   +           174                         0                                  3   …  
#   Interval_4     1   3364506   3364805   +           300                         5                                  1   …  
#   Interval_5     1   3445884   3446410   +           527                         2                                  5   …  
#   Interval_6     1   3475110   3475305   +           196                         2                                  2   …  
#   Interval_7     1   3506104   3506325   +           222                         0                                  3   …  
#   Interval_8     1   3564768   3564984   +           217                         0                                  5   …  
#  ------------ ----- --------- --------- -------- -------- ------------------------- ---------------------------------- --- 


# 1.2) Boolean consensus peaks matrix
# +---------------------------------+------------+------------+
# | chr                             | 1          | 1          |
# | start                           | 3008719    | 3009020    |
# | end                             | 3042821    | 3043241    |
# | interval_id                     | Interval_1 | Interval_2 |
# | num_peaks                       | 1          | 6          |
# | num_samples                     | 1          | 6          |
# | BcellBoneMarrow.mLb.clN.bool    | FALSE      | FALSE      |
# | CD4SinglePositive.mLb.clN.bool  | FALSE      | FALSE      |
# | …                               | …          | …          |
# | BcellBoneMarrow.mLb.clN.fc      | NA         | NA         |
# | CD4SinglePositive.mLb.clN.fc    | NA         | NA         |
# | …                               | …          | …          |
# | BcellBoneMarrow.mLb.clN.qval    | NA         | NA         |
# | CD4SinglePositive.mLb.clN.qval  | NA         | NA         |
# | …                               | …          | …          |
# | BcellBoneMarrow.mLb.clN.pval    | NA         | NA         |
# | CD4SinglePositive.mLb.clN.pval  | NA         | NA         |
# | …                               | …          | …          |
# | BcellBoneMarrow.mLb.clN.start   | NA         | NA         |
# | CD4SinglePositive.mLb.clN.start | NA         | NA         |
# | …                               | …          | …          |
# | BcellBoneMarrow.mLb.clN.end     | NA         | NA         |
# | CD4SinglePositive.mLb.clN.end   | NA         | NA         |
# +---------------------------------+------------+------------+

