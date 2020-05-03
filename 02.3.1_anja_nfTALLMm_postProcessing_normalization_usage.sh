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
suppressPackageStartupMessages(library("ggplot2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("reshape2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("pander", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("Hmisc", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts=FALSE, quietly=TRUE))

inputFile  <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_analysis_consensus_peaks_rawCounts.txt"
annotFile  <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_sample_annotation.txt"
outputFile <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_analysis_consensus_peaks_vstCounts.txt"

# Import data from featureCounts
countdata   <- fread(inputFile, header=TRUE)
sampleTable <- fread(annotFile, header=TRUE)

# Filter out peaks where not a single sample has more than 50 reads.
countdata = countdata[apply(countdata[,5:ncol(countdata)], 1, max) > 50,]

# Subset genomic ranges from the counts DT
grangesDT  <- countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID')]
# Filter out peaks where not a single sample has more than 50 reads.
grangesDT <- grangesDT[apply(countdata[,4:ncol(countdata)], 1, max) > 50,]

# Remove PeakChrom PeakStart  PeakEnd columns
countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID'):=NULL]

# Get the dds object
# the design formula contains a numeric variable with integer values,
# specifying a model with increasing fold change for higher values.
# did you mean for this to be a factor? if so, first convert
# this variable to a factor using the factor() function
sampleTable$Study <- as.factor(sampleTable$Study)
dds <- DESeqDataSetFromMatrix(colData  = sampleTable, countData=countdata, design =~Study)

# VST 
vst      <- varianceStabilizingTransformation(dds, blind=FALSE)  
vsd      <- assay(vst)

# Add the genomic ranges in base R
vstBindDT <- cbind(grangesDT, vsd)

# Round all normalized counts to 4 decimal places
vstOutDT  <- data.table(vstBindDT %>% mutate_at(vars(-PeakChrom,-PeakStart,-PeakEnd,-PeakID), funs(round(., 4))))

# Save the counts file to output csv file
fwrite(vstOutDT   , file = outputFile , sep = '\t', quote = F)

# PCA
pca = prcomp(t(vsd))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
si <- as.data.frame(sampleTable)
rownames(si) <- si$SampleName; si$SampleName <- NULL
pcaData=merge(pcaData, si, by.x=0, by.y=0)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, group=CellLine)) + geom_point(aes(size=3, shape=Study, color=CellLine))
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
normPlotPDF  <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_vst_norm_PCA.pdf"
ggsave(normPlotPDF, width=15)

# Data visualization of normalized counts
lf = melt(vsd, id.vars=c())
pander(head(lf))
pg1 <- ggplot(data=lf, aes(x=Var2, y=value)) + geom_boxplot(aes(group=Var2)) + xlab("SampleName") + ylab("Normalized Count") + coord_flip()
normPlotPDF  <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_vst_sample_boxplot.pdf"
ggsave(normPlotPDF)
pg2 <- ggplot(data=lf, aes(x=value, after_stat(density))) + geom_freqpoly(aes(group=Var2, color=Var2), bins=30) + xlab("Normalized Count") + ylab("density")
normPlotPDF  <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_vst_normalization.pdf"
ggsave(normPlotPDF, width=15)

# Sample cluster dendrogram
hclustPlotPDF  <- "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_vst_sample_hclust_dendrogram.pdf"
dists          <- dist(t(assay(vst)))
pdf(hclustPlotPDF); par(mar=c(20,4,1,1));
plot(hclust(dists))
dev.off()

# Heatmap of sample-to-sample distances using the Poisson Distance
suppressMessages(library("PoiClaClu"))
suppressMessages(library("pheatmap"))
poisd <- PoissonDistance(t(vsd))

# In order to plot the sample distance matrix with the rows/columns arranged by the distances in our distance matrix, we manually provide sampleDists 
# to the clustering_distance argument of the pheatmap function. Otherwise the pheatmap function would assume that the matrix contains the data values 
# themselves, and would calculate distances between the rows/columns of the distance matrix, which is not desired.
samplePoisDistMatrix <- as.matrix( poisd$dd)

# Change the row names of the distance matrix to contain Treatment and Control instead of sample ID, so that we have all this information in view when looking at the heatmap.
rownames(samplePoisDistMatrix) <- paste( attr(rld,"colData")$condition, rownames(attr(rld,"colData")),sep="-")
colnames(samplePoisDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
# colors <- colorRampPalette(rev(c('gold','darkorange','darkred')))(256)
# colors <- colorRampPalette(c('darkblue','blue','lightblue','white','orange','red','darkred'))(1024)
pheatmap(samplePoisDistMatrix,
        clustering_distance_rows=poisd$dd,
        clustering_distance_cols=poisd$dd,
        col=colors,
        fontsize=8)

#par(new=TRUE)
# Turn off device driver (to flush output to PNG/PDF file)
dev.off()
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

