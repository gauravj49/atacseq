# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Source: https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
# For later purpose
# https://bioconductor.org/packages/release/bioc/vignettes/rGREAT/inst/doc/rGREAT.html

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

# Annotate each sample peaks with ChIPseeker
library(genomation)
library(ChIPseeker)
library(UpSetR)
library(ReactomePA)
library(clusterProfiler)
library(rGREAT)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene 

# https://rdrr.io/bioc/GenomicFeatures/man/makeTxDbFromGFF.html

outdir  <- '/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/peaks_annotation'
system(paste("mkdir -p ", outdir, sep=''))

# List of sample peaks
peaksFileList <- system("ls /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/*/macs2peaks/*.bed", intern = TRUE)

# 1) ============================================
# 1) Annotate and visualize peaks
# Loop doesn`t work here so doing it one by one
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup/macs2peaks/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup_summits.bed"                
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_LivMet-3_005_atac_031_S12_R1_001_rmdup/macs2peaks/5320_LivMet-3_005_atac_031_S12_R1_001_rmdup_summits.bed"                
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_LungMet-1_005_atac_029_S10_R1_001_rmdup/macs2peaks/5320_LungMet-1_005_atac_029_S10_R1_001_rmdup_summits.bed"              
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_PPT-1_004_atac_023_20k_S18_R1_001_rmdup/macs2peaks/5320_PPT-1_004_atac_023_20k_S18_R1_001_rmdup_summits.bed"              
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_PPT-1_005_atac_028_S9_R1_001_rmdup/macs2peaks/5320_PPT-1_005_atac_028_S9_R1_001_rmdup_summits.bed"                        
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_LivMet-1_005_atac_025_S6_R1_001_rmdup/macs2peaks/53631_LivMet-1_005_atac_025_S6_R1_001_rmdup_summits.bed"                
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_LungMet-2_005_atac_026_S7_R1_001_rmdup/macs2peaks/53631_LungMet-2_005_atac_026_S7_R1_001_rmdup_summits.bed"              
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_LungMet-3_005_atac_027_S8_R1_001_rmdup/macs2peaks/53631_LungMet-3_005_atac_027_S8_R1_001_rmdup_summits.bed"              
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_PPT-1_005_atac_024_S5_R1_001_rmdup/macs2peaks/53631_PPT-1_005_atac_024_S5_R1_001_rmdup_summits.bed"                      
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_LivMet-1_005_atac_033_S2_R1_001_rmdup/macs2peaks/53646_LivMet-1_005_atac_033_S2_R1_001_rmdup_summits.bed"                
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_LivMet-2_005_atac_034_S3_R1_001_rmdup/macs2peaks/53646_LivMet-2_005_atac_034_S3_R1_001_rmdup_summits.bed"                
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_LivMet-3_005_atac_035_S4_R1_001_rmdup/macs2peaks/53646_LivMet-3_005_atac_035_S4_R1_001_rmdup_summits.bed"                
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_PPT-1_004_atac_022_S17_R1_001_rmdup/macs2peaks/53646_PPT-1_004_atac_022_S17_R1_001_rmdup_summits.bed"                    
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_PPT-1_005_atac_032_S1_R1_001_rmdup/macs2peaks/53646_PPT-1_005_atac_032_S1_R1_001_rmdup_summits.bed"                      
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_PPT-1_frozen_004_atac_022_fr_S16_R1_001_rmdup/macs2peaks/53646_PPT-1_frozen_004_atac_022_fr_S16_R1_001_rmdup_summits.bed"
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/6075_blood_005_atac_038_S15_R1_001_rmdup/macs2peaks/6075_blood_005_atac_038_S15_R1_001_rmdup_summits.bed"                      
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/6075_LungMet-1_005_atac_037_S14_R1_001_rmdup/macs2peaks/6075_LungMet-1_005_atac_037_S14_R1_001_rmdup_summits.bed"              
# p <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/6075_PPT-1_005_atac_036_S13_R1_001_rmdup/macs2peaks/6075_PPT-1_005_atac_036_S13_R1_001_rmdup_summits.bed" 

indPeaksAnnDir <- paste0(outdir,'/individual_samples'); system(paste("mkdir -p ", indPeaksAnnDir, sep=''))
# Get the base file name
bname <- basename(tools::file_path_sans_ext(p))
cat("\n- Processing sample: ", bname, "\n")

# Read bed file
cat("\n\t- Reading the bed file\n")
peaksBed <- readBed(p, track.line=FALSE,remove.unusual=TRUE)

# Annotate regions
cat("\n\t- Annotate regions\n")
peakAnno <- annotatePeak(peaksBed, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# Visualize genomic annotation
cat("\n\t- Visualize genomic annotation\n")
peaksAnnPdf  <- paste0(indPeaksAnnDir, '/',bname,'_chipseeker_annotation_plots.pdf', sep='')
# get_ann_plots(peakAnno, peaksAnnPdf)
pdf(peaksAnnPdf)
cat("\nPlotting AnnoPie")
plotAnnoPie(peakAnno)
cat("\nPlotting plotAnnoBar")
plotAnnoBar(peakAnno)
cat("\nPlotting upsetplot")
upsetplot(peakAnno)
dev.off()

# 2) ============================================
# 2) ChIP peak data set comparison
# # 2.1) Profile of several ChIP peak data binding to TSS region
# promoter      <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
# tagMatrixList <- lapply(peaksFileList, getTagMatrix, windows=promoter)
# tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL) # Gives an error (Error in plot.new() : figure margins too large)

# 3) ============================================
# 3) ChIP peak annotation comparision
peakAnnoList <- lapply(peaksFileList, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

# Add the names to the list
names(peakAnnoList) <- basename(tools::file_path_sans_ext(peaksFileList))
compPeaksAnnDir <- paste0(outdir,'/sampleComparisons'); system(paste("mkdir -p ", compPeaksAnnDir, sep=''))

# Genomic annotation among different samples
cat("\n\t- Visualize genomic annotation among different samples\n")
peaksAnnPdf  <- paste0(compPeaksAnnDir, '/',bname,'_chipseeker_genomicAnnotation_all_samples.pdf', sep='')
pdf(peaksAnnPdf, width = 25, height = 10); plotAnnoBar(peakAnnoList); dev.off()

# Distribution of binding sites among different samples
cat("\n\t- Visualize distribution of binding sites among different samples\n")
peaksAnnPdf  <- paste0(compPeaksAnnDir, '/',bname,'_chipseeker_bindingSitesDistribution_all_samples.pdf', sep='')
pdf(peaksAnnPdf, width = 25, height = 10); plotDistToTSS(peakAnnoList); dev.off()

# 4) ============================================
# 4) Functional profiles comparison
genes        <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# genes <- lapply(genes, unique) 
names(genes) <- sub("_", "\n", names(genes))
compKEGG     <- compareCluster(geneCluster   = genes,
                               fun           = "enrichKEGG",
                               pvalueCutoff  = 0.5,
                               pAdjustMethod = "BH")
KPEAPdf  <- paste0(compPeaksAnnDir, '/',bname,'_chipseeker_KEGGPathwayEnrichmentAnalysis_all_samples.pdf', sep='')
pdf(peaksAnnPdf, width = 25, height = 10);
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis"); dev.off()

# 5) ============================================
# 5) Functional profiles comparison with GREAT
#    https://github.com/jokergoo/rGREAT
set.seed(123)

job = submitGreatJob(peaksBed, species="mm10")
job

tb = getEnrichmentTables(job)
names(tb)


# Finding Enriched Motifs in Genomic Regions (findMotifsGenome.pl)
inputDir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/homer/input/Top1pcATACpeaks_Homer"
species="mm10"
for f in ATAC001_heatmap_cluster1.bed  ATAC001_heatmap_cluster3.bed  ATAC001_heatmap_cluster5.bed  ATAC001_heatmap_cluster7.bed ATAC001_heatmap_cluster2.bed  ATAC001_heatmap_cluster4.bed  ATAC001_heatmap_cluster6.bed  ATAC001_heatmap_cluster8.bed;
do
    motifoutdir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/homer/output/motifs/$(basename ${f} .bed)"; mkdir -p ${motifoutdir}
    echo "findMotifsGenome.pl ${inputDir}/${f} ${species} ${motifoutdir} -len 6,8,10 -p 8"
    findMotifsGenome.pl ${inputDir}/${f} ${species} ${motifoutdir} -len 6,8,10 -p 8
done;

