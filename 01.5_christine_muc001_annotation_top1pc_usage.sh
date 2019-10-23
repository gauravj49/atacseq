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

# Annotate each sample peaks with ChIPseeker
library(genomation)
library(ChIPseeker)
library(UpSetR)
library(ReactomePA)
library(clusterProfiler)
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


for (p in peaksFileList){
    # Get the base file name
    bname <- basename(tools::file_path_sans_ext(p))
    cat("\n- Processing sample: ", bname, "\n")

    # Read bed file
    cat("\n\t- Reading the bed file\n")
    peaksBed <- readBed(p, track.line=FALSE,remove.unusual=TRUE)

    # Annotate regions
    cat("\n\t- Annotate regions\n")
    peakAnno <- annotatePeak(peaksBed, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

    pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
    head(pathway1, 2)

}
############ USER DEFINED FUNCTIONS ##########
# Generate RLE plots 
get_ann_plots <- function(peakAnno, peaksAnnPdf){
    pdf(peaksAnnPdf)
    cat("\nPlotting AnnoPie")
    plotAnnoPie(peakAnno)
    cat("\nPlotting plotAnnoBar")
    plotAnnoBar(peakAnno)
    cat("\nPlotting upsetplot")
    upsetplot(peakAnno)
}


# https://rdrr.io/bioc/genomation/man/readBed.html
allpirs <- readBed('/fshare/users/jaing/bin/projects/DeAnalysis/input/annotation/bed/hg19_pirnaBank_piRs.bed',track.line=FALSE,remove.unusual=TRUE)
csfpirs <- readBed('/fshare/users/jaing/bin/projects/DeAnalysis/input/annotation/bed/hg19_pirnaBank_merged_cohort12_all_samples_ab450tau200_diagnosis_filtered_pirnaome.bed',track.line=FALSE,remove.unusual=TRUE)
sigpirs <- readBed('/fshare/users/jaing/bin/projects/DeAnalysis/input/annotation/bed/signature_hg19_pirnaBank_merged_cohort12_all_samples_ab450tau200_diagnosis_filtered_pirnaome.bed',track.line=FALSE,remove.unusual=TRUE)

# Annotate regions
# 1) https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
# 2) https://www.rdocumentation.org/packages/ChIPseeker/versions/1.8.6/topics/annotatePeak
allpirAnno <- annotatePeak(allpirs, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
csfpirAnno <- annotatePeak(csfpirs, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
sigpirAnno <- annotatePeak(sigpirs, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

#### For all piRNAome ####
# Visualize Genomic Annotation
allPirPdf  <- paste(outdir, '/hg19_pirnaBank_piRs_annotation_plots.pdf', sep='')
pdf(allPirPdf)
plotAnnoPie(allpirAnno)
plotAnnoBar(allpirAnno)
vennpie(allpirAnno)
dev.off()
# We can combine vennpie with upsetplot by setting vennpie = TRUE.
allPirPng  <- paste(outdir, '/hg19_pirnaBank_piRs_annotation_upsetplot.png', sep='')
png(filename=allPirPng, height=3000, width=6000, res=300)
upsetplot(allpirAnno, vennpie=TRUE)
dev.off()

#### For csf piRNAs ####
# Visualize Genomic Annotation
csfPirPdf  <- paste(outdir, '/hg19_pirnaBank_merged_cohort12_filtered_pimirnaome_plots.pdf', sep='')
pdf(csfPirPdf)
plotAnnoPie(csfpirAnno)
plotAnnoBar(csfpirAnno)
vennpie(csfpirAnno)
dev.off()
# We can combine vennpie with upsetplot by setting vennpie = TRUE.
csfPirPng  <- paste(outdir, '/hg19_pirnaBank_merged_cohort12_filtered_pimirnaome_upsetplot.png', sep='')
png(filename=csfPirPng, height=3000, width=6000, res=300)
upsetplot(csfpirAnno, vennpie=TRUE)
dev.off()


#### For csignaturesf piRNAs ####
# Visualize Genomic Annotation
sigPirPdf  <- paste(outdir, '/signature_hg19_pirnaBank_merged_cohort12_filtered_pimirnaome_plots.pdf', sep='')
pdf(sigPirPdf)
plotAnnoPie(sigpirAnno)
plotAnnoBar(sigpirAnno)
vennpie(sigpirAnno)
dev.off()
# We can combine vennpie with upsetplot by setting vennpie = TRUE.
sigPirPng  <- paste(outdir, '/signature_hg19_pirnaBank_merged_cohort12_filtered_pimirnaome_upsetplot.png', sep='')
png(filename=sigPirPng, height=3000, width=6000, res=300)
upsetplot(sigpirAnno, vennpie=TRUE)
dev.off()

# Save the signature annotation in the tab separated text file
sigPirtxt  <- paste(outdir, '/signature_hg19_pirnaBank_merged_cohort12_filtered_pimirnaome.txt', sep='')
# Commented and uncommented statement provides the same result
# write.table(sigpirAnno@anno@elementMetadata@listData, file=sigPirtxt, quote=FALSE, sep="\t")
write.table(sigpirAnno, file=sigPirtxt, quote=FALSE, sep="\t")
















































