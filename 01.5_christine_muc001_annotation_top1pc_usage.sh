# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Source: https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
# For later purpose
https://bioconductor.org/packages/release/bioc/vignettes/rGREAT/inst/doc/rGREAT.html

R 
# Install R packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ChIPseeker")
# BiocManager::install("genomation")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("UpSetR")

# Annotate each sample peaks with ChIPseeker
library(genomation)
library(ChIPseeker)
library(UpSetR)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

outdir  <- '/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/peaks_annotation'
system(paste("mkdir -p ", outdir, sep=''))

# List of sample peaks
peaksFileList <- system("ls /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/*/macs2peaks/*.bed", intern = TRUE)

# Annotate and visualize peaks
indPeaksAnnDir <- paste0(outdir,'/individual_samples'); system(paste("mkdir -p ", indPeaksAnnDir, sep=''))
for (p in peaksFileList){
    # Get the base file name
    bname <- basename(tools::file_path_sans_ext(p))
    cat("\n- Processing sample: ", bname, "\n")

    # Read bed file
    cat("\n\t- Reading the bed file\n")
    peaksBed <- readBed(p, track.line=FALSE,remove.unusual=TRUE)

    # Annotate regions
    cat("\n\t- Annotate regions\n")
    peaksAnn <- annotatePeak(peaksBed, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

    # Visualize genomic annotation
    cat("\n\t- Visualize genomic annotation\n")
    peaksAnnPdf  <- paste0(indPeaksAnnDir, '/',bname,'_chipseeker_annotation_plots.pdf', sep='')
    pdf(peaksAnnPdf)
    plotAnnoPie(peaksAnn)
    plotAnnoBar(peaksAnn)
    upsetplot(peakAnno)
    dev.off()
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
