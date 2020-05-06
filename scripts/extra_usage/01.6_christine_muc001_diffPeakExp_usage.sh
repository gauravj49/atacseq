# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Source: https://tobiasrausch.com/courses/atac/atac-seq-data-analysis.html#r-libraries
R 
# Install R packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ChIPseeker")
# BiocManager::install("genomation")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("UpSetR")
# BiocManager::install("ReactomePA")
# BiocManager::install("clusterProfiler")
# BiocManager::install("rGREAT")
# BiocManager::install("pander")
# BiocManager::install("pastecs")

# Annotate each sample peaks with ChIPseeker
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(rGREAT))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.ensGene))
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene 

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(pastecs))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(pheatmap))

# 1) Load merged ATAC-Seq count matrix and the sample information 
inputPeaksFile  <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.matrix"
sample_info     <- "/home/rad/users/gaurav/projects/seqAnalysis/atacseq/docs/samples_info.txt"
outdir          <- '/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/peaks_annotation'
compPeaksAnnDir <- paste0(outdir,'/sampleComparisons'); system(paste("mkdir -p ", compPeaksAnnDir, sep=''))

# Get the base file name
bname <- basename(tools::file_path_sans_ext(inputPeaksFile)); bname <-gsub("_rawCounts", "", x)
cat("\n- Processing sample: ", bname, "\n")

# Import data from featureCounts
countdataDF  <- read.table(inputPeaksFile, header=TRUE, row.names=1, check.names = FALSE)
# ┌─────────────────────────────────────┬────────────────────────────────┬────────────────────────────────┬─────┐
# │               peakID                │ 5320_livmet-1_005_atac_030_s11 │ 5320_livmet-3_005_atac_031_s12 │ ... │
# ├─────────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼─────┤
# │ chr2_104603717_104604218_atacPeak_0 │                              7 │                              7 │ ... │
# │ chr2_104607403_104607904_atacPeak_1 │                              6 │                             21 │ ... │
# │ ...                                 │                            ... │                            ... │     │
# └─────────────────────────────────────┴────────────────────────────────┴────────────────────────────────┴─────┘

sampleInfoDF         <- read.table(sample_info, header=TRUE, sep='\t', row.names=1)
sampleInfoDF$donor   <- factor(sampleInfoDF$donor)
sampleInfoDF$cluster <- factor(sampleInfoDF$cluster)
sampleInfoDF$Organ   <- factor(sampleInfoDF$organ)
# ┌───────────────────────────────────────┬─────────┬───────┬────────┐
# │                sample                 │ cluster │ donor │ Organ  │
# ├───────────────────────────────────────┼─────────┼───────┼────────┤
# │ 5320_LivMet-1_005_atac_030_S11_R1_001 │ C1      │  5320 │ LivMet │
# │ 5320_LivMet-3_005_atac_031_S12_R1_001 │ C1      │  5320 │ LivMet │
# │ ...                                   │ ...     │   ... │ ...    │
# └───────────────────────────────────────┴─────────┴───────┴────────┘



# df = read.table("atac.data.gz", header=T)
# si = read.table("blood.samples", header=F)
# colnames(si) = c("sample", "celltype", "donor")
# rownames(si) = si$sample
# si$donor = factor(si$donor)
pander(dim(countdataDF), "Data dimensions")  # _159990_ and _18_
print(summary(sampleInfoDF))
#  cluster    donor       Organ  
#  C1 : 5   5320 :5   blood  :1  
#  C2b:10   6075 :3   LivMet :6  
#  C2c: 3   53631:4   LungMet:4  
#           53646:6   PPT    :7

# 2) Removing missing peaks
countdataDF <- countdataDF[apply(countdataDF[,4:ncol(countdataDF)], 1, max) > 50,]
pander(dim(countdataDF), "Data dimensions") # # _90724_ and _18_

# 2.2) Data Exploration
# The other major problem in genomic count data sets is that they often show heteroscedasticity which means in our case that different peaks show different levels of variabilities in the number of reads. This is a major problem for differential peak calling.
rowsummary = data.frame(rowmeans = apply(countdataDF[, 4:ncol(countdataDF)], 1, mean), rowsds = apply(countdataDF[, 4:ncol(countdataDF)], 1, sd))

sdMeanPdf  <- paste0(compPeaksAnnDir, '/',bname,'_SD_Mean_all_samples.png', sep='')
png(sdMeanPdf); 
ggplot(data=rowsummary, aes(x=rowmeans, y=rowsds)) + geom_point() + xlab("Peak means") + ylab("Peak SDs")
dev.off()

# 3) Differential peak calling 
This line is need to fix Error in .validate_names(colnames, ans_colnames, "assay colnames()", "colData rownames()") : 
colnames(countdataDF) <- NULL 

dds <- DESeqDataSetFromMatrix(countData = countdataDF, colData = sampleInfoDF, design = ~donor)
dds <- DESeq(dds, fitType='local')

resultsNames(dds)

resdonor6075vs5320  <- results(dds, lfcThreshold=1.5, contrast=c("donor", "6075", "5320"))
resdonor53631vs5320 <- results(dds, lfcThreshold=1.5, contrast=c("donor", "53631", "5320"))
resdonor53646vs5320 <- results(dds, lfcThreshold=1.5, contrast=c("donor", "53646", "5320"))

# print(mcols(res, use.names=T))
print(summary(resdonor6075vs5320))
print(summary(resdonor53631vs5320))
print(summary(resdonor53646vs5320))



hist(res$pvalue, breaks=0:20/20, col="grey50", border="white", xlim=c(0,1), main="Histogram of p-values", xlab="p-value")

# Data normalization
vst    <- varianceStabilizingTransformation(dds, fitType='local')  
vsd    <- assay(vst)
normDF <- data.frame(vsd, check.names = FALSE)

# print(sum(abs(res$log2FoldChange) > 2))
# Heatmaps
mat = normDF[which(abs(resdonor6075vs5320$log2FoldChange) >2),]
mat = mat - rowMeans(mat)
anno = as.data.frame(colData(dds)[, c("donor", "cluster")])
rownames(mat) = NULL
heatmapPdf  <- paste0(compPeaksAnnDir, '/',bname,'_DEPeaks_Heatmap_6075vs5320.pdf', sep='')
pheatmap(mat, annotation_col = anno, scale="row", filename = heatmapPdf, show_rownames=FALSE)

mat = normDF[which(abs(resdonor53631vs5320$log2FoldChange) >2),]
mat = mat - rowMeans(mat)
anno = as.data.frame(colData(dds)[, c("donor", "cluster")])
rownames(mat) = NULL
heatmapPdf  <- paste0(compPeaksAnnDir, '/',bname,'_DEPeaks_Heatmap_53631vs5320.pdf', sep='')
pheatmap(mat, annotation_col = anno, scale="row", filename = heatmapPdf, show_rownames=FALSE)

mat = normDF[which(abs(resdonor53646vs5320$log2FoldChange) >2),]
mat = mat - rowMeans(mat)
anno = as.data.frame(colData(dds)[, c("donor", "cluster")])
rownames(mat) = NULL
heatmapPdf  <- paste0(compPeaksAnnDir, '/',bname,'_DEPeaks_Heatmap_53646vs5320.pdf', sep='')
pheatmap(mat, annotation_col = anno, scale="row", filename = heatmapPdf, show_rownames=FALSE)

# 4) Data visualization
lf = melt(normDF, id.vars=c())
pander(head(lf))
vstPdf  <- paste0(compPeaksAnnDir, '/',bname,'_VST_all_samples.pdf', sep='')
pdf(vstPdf, width = 8, height = 10); 
ggplot(data=lf, aes(x=Var2, y=value)) + geom_boxplot(aes(group=Var2)) + xlab("Sample IDs") + ylab("VST Normalized Count") + coord_flip()
dev.off()


# 4.2) PCA
pca = prcomp(t(normDF))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData=cbind(pcaData, sampleInfoDF)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
pcaPdf  <- paste0(compPeaksAnnDir, '/',bname,'_PCA_all_samples.pdf', sep='')
pdf(pcaPdf, width = 10, height = 8); 
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, color=donor, shape=Organ, size=cluster)) + geom_point(size=3)
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
print(p)
dev.off()

# 4.3) PCA Loadings
# Inspect the loadings for each PC to know which peaks contribute most to the separation of the individual clusters
loadings = abs(pca$rotation)
contribution = as.data.frame(sweep(loadings, 2, colSums(loadings), "/"))
contribution = contribution[with(contribution, order(-PC1)),]
pander(head(contribution))



