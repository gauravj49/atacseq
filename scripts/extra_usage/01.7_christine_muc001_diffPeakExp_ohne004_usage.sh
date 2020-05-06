# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

R 
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
outdir          <- '/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/peaks_annotation/ohne004'
compPeaksAnnDir <- paste0(outdir,'/sampleComparisons'); system(paste("mkdir -p ", compPeaksAnnDir, sep=''))

# Get the base file name
bname <- basename(tools::file_path_sans_ext(inputPeaksFile)); bname <-gsub("_rawCounts", "", x)
cat("\n- Processing sample: ", bname, "\n")

# Import data from featureCounts
countdataDF  <- read.table(inputPeaksFile, header=TRUE, row.names=1, check.names = FALSE)
countdataDF  <- countdataDF[, -grep("004", colnames(countdataDF))]
# ┌─────────────────────────────────────┬────────────────────────────────┬────────────────────────────────┬─────┐
# │               peakID                │ 5320_livmet-1_005_atac_030_s11 │ 5320_livmet-3_005_atac_031_s12 │ ... │
# ├─────────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼─────┤
# │ chr2_104603717_104604218_atacPeak_0 │                              7 │                              7 │ ... │
# │ chr2_104607403_104607904_atacPeak_1 │                              6 │                             21 │ ... │
# │ ...                                 │                            ... │                            ... │     │
# └─────────────────────────────────────┴────────────────────────────────┴────────────────────────────────┴─────┘

sampleInfoDF         <- read.table(sample_info, header=TRUE, sep='\t', row.names=1)
sampleInfoDF         <- sampleInfoDF[-grep("004",rownames(sampleInfoDF)),]
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
pander(dim(countdataDF), "Data dimensions")  # _159990_ and _15_
print(summary(sampleInfoDF))
#  cluster   donor       Organ  
#  C1 :4   5320 :4   blood  :1  
#  C2b:8   6075 :3   LivMet :6  
#  C2c:3   53631:4   LungMet:4  
#          53646:4   PPT    :4 

# 2) Removing missing peaks
countdataDF <- countdataDF[apply(countdataDF[,4:ncol(countdataDF)], 1, max) > 50,]
pander(dim(countdataDF), "Data dimensions") # # _90724_ and _18_

# 2.2) Data Exploration
# The other major problem in genomic count data sets is that they often show heteroscedasticity which means in our case that different peaks show different levels of variabilities in the number of reads. This is a major problem for differential peak calling.
rowsummary = data.frame(rowmeans = apply(countdataDF[, 4:ncol(countdataDF)], 1, mean), rowsds = apply(countdataDF[, 4:ncol(countdataDF)], 1, sd))

sdMeanPdf  <- paste0(compPeaksAnnDir, '/',bname,'_SD_Mean_ohne004_samples.png', sep='')
png(sdMeanPdf); 
ggplot(data=rowsummary, aes(x=rowmeans, y=rowsds)) + geom_point() + xlab("Peak means") + ylab("Peak SDs")
dev.off()

# 3) Differential peak calling 
# This line is need to fix Error in .validate_names(colnames, ans_colnames, "assay colnames()", "colData rownames()") : 
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

resdonor6075vs5320  <- resdonor6075vs5320[order(resdonor6075vs5320$padj),   ]
resdonor53631vs5320 <- resdonor53631vs5320[order(resdonor53631vs5320$padj), ]
resdonor53646vs5320 <- resdonor53646vs5320[order(resdonor53646vs5320$padj), ]

# Data normalization
vst    <- varianceStabilizingTransformation(dds, fitType='local')  
vsd    <- assay(vst)
normDF <- data.frame(vsd, check.names = FALSE)

# Save data results and normalized reads to csv!
resultsDir      <- paste(outdir, "/diffPeaks", sep=''); system(paste("mkdir -p ", resultsDir, sep=''))
cat("\t7.2) Save data results and normalized reads to csv ...\n")
# resdata <- merge(as.data.frame(res), as.data.frame(counts(fdds,normalized=T)), by='row.names',sort=F)
resultsOutFile  <- paste(resultsDir, "/", bname,"_DEPeaks_resdonor6075vs5320_RESULTS.txt", sep='')
resdata <- merge(as.data.frame(resdonor6075vs5320), normDF, by='row.names',sort=F)
names(resdata)[1] <- 'feature'
write.table(resdata, file = resultsOutFile, row.names = F, sep = '\t', quote = F)

resultsOutFile  <- paste(resultsDir, "/", bname,"_DEPeaks_resdonor53631vs5320_RESULTS.txt", sep='')
resdata <- merge(as.data.frame(resdonor53631vs5320), normDF, by='row.names',sort=F)
names(resdata)[1] <- 'feature'
write.table(resdata, file = resultsOutFile, row.names = F, sep = '\t', quote = F)

resultsOutFile  <- paste(resultsDir, "/", bname,"_DEPeaks_resdonor53646vs5320_RESULTS.txt", sep='')
resdata <- merge(as.data.frame(resdonor53646vs5320), normDF, by='row.names',sort=F)
names(resdata)[1] <- 'feature'
write.table(resdata, file = resultsOutFile, row.names = F, sep = '\t', quote = F)


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
vstPdf  <- paste0(compPeaksAnnDir, '/',bname,'_VST_ohne004_samples.pdf', sep='')
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
pcaPdf  <- paste0(compPeaksAnnDir, '/',bname,'_PCA_ohne004_samples.pdf', sep='')
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



suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts=FALSE, quietly=TRUE))
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(factor(sampleCondition)))])
cat("\n8) Plotting diagnosis plots\n")
cat("\t8.1) Plotting PCA plots...\n")
outputFile <- paste0(outdir,'/sampleComparisons/',bname,'_PCA_ohne004_samples_final')
plotpca(vsd, compPeaksAnnDir, outputFile, mycols)

cat("\t8.2) Plotting volcano plots...\n")
tryCatch({
plot_volcanoPlot(resdata, digPlotDirPng, digPlotDirPdf, outputFile)
},
error=function(e){
	cat("\n\t- Error in plotting volcano plots")
	print(e)
},
warning=function(w){
	print(w)
})

cat("\t8.3) Plotting MA plots...\n")
tryCatch({
plot_maPlot(resdata, digPlotDirPng, digPlotDirPdf, outputFile)
},
error=function(e){
	cat("\n\t- Error in plotting MA plots")
	print(e)
},
warning=function(w){
	print(w)
})

cat("\t8.4) Plotting pvalues plots...\n")
tryCatch({
plot_pvalPlot(res, digPlotDirPng, digPlotDirPdf, outputFile)
},
error=function(e){
	cat("\n\t- Error in plotting pvalue plots:")
	print(e)
},
warning=function(w){
	print(w)
})

cat("\t8.5) Plotting independent filtering plots...\n")
tryCatch({
	plot_indfilPlot(res, digPlotDirPng, digPlotDirPdf, outputFile)
},
error=function(e){
	cat("\n\t- Error in plotting independent filtering plots")
	print(e)
},
warning=function(w){
	print(w)
})

tryCatch({
	cat("\t8.6) Plotting sample correlation heatmaps...\n")
plot_sampleDistplot(fdds, rld, digPlotDirPng, digPlotDirPdf, outputFile)
},
error=function(e){
	cat("\n\t- Error in plotting sample correlation heatmaps")
	print(e)
},
warning=function(w){
	print(w)
})


cat("\n ***********************************\n")
cat("- Output Results Files          : ", resultsOutFile, "\n")
options(warn=0)

########## USER DEFINED FUNCTIONS #################
# Draw PCA plots
plotpca <- function(rld, pcaDir, outputFile, mycols){
	# Create the png filename
	pngfile  <- paste(pcaDir, "/", basename(outputFile), "_pca.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	rld_pca(rld, colors=mycols, intgroup="x", main=paste("PCA (", conditionVals[2],", ", conditionVals[1],")", sep=""), pcaDir=pcaDir)

	# Create the pdf filename
	pdffile<-paste(pcaDir, "/", basename(outputFile), "_pca.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	rld_pca(rld, colors=mycols, intgroup="x", main=paste("PCA (", conditionVals[2],", ", conditionVals[1],")", sep=""), plotlabels=TRUE, pcaDir=pcaDir)
}

rld_pca <- function (rld, intgroup = "condition", ntop = 500, plotlabels=FALSE, colors=NULL, legendpos="topright", main="PCA Biplot", textcx=1, pcaDir, ...) {
	suppressPackageStartupMessages(require(genefilter, warn.conflicts=FALSE, quietly=TRUE))
	suppressPackageStartupMessages(require(calibrate, warn.conflicts=FALSE, quietly=TRUE))
	suppressPackageStartupMessages(require(RColorBrewer, warn.conflicts=FALSE, quietly=TRUE))
	rv = rowVars(assay(rld))
	select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	pca = prcomp(t(assay(rld)[select, ]))
	fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
	if (is.null(colors)) {
		if (nlevels(fac) >= 3) {
			colors = brewer.pal(nlevels(fac), "Paired")
		} else {
			colors = c("black", "red")
		}
	}
	
	cexsize = 2.0
	if (plotlabels){
		# Get the names
		rnames          <- rownames(as.data.frame(pca$x))
		nnum            <- seq(1:length(rnames))
		rarray          <- setNames(nnum, rnames)
		pcIDs           <- as.data.frame(as.matrix(rarray))

		# Rename first column as snow
		names(pcIDs)    <- c("sno")
		
		# Put rownames as second column
 		pcIDs <- cbind(data.frame(pcIDs, row.names=NULL), rownames(pcIDs))

 		# Rename the second column
		colnames(pcIDs)[2] <- 'sampleNames'

		# Save the dataframe in the output file
		txtFiledir      <- paste(pcaDir,"/fileIds", sep=''); system(paste("mkdir -p ", txtFiledir, sep=''));
		sampleIDfile    <- paste(txtFiledir, "/", basename(outputFile), "_fileIds.txt",sep='');
		write.table(pcIDs, file = sampleIDfile, row.names = F, sep = '\t', quote = F)
		
	
		# Put legend on bottom 1/8th of the chart
		layout(matrix(c(1,2,2)), heights=c(length(rnames)*6, length(rnames)*2))  

		# circle size
		cexsize = 1.0 
	}

	# Plot the pca
	pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
	pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
	pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
	pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
	
	plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], col=NULL, pch=21, cex=cexsize, xlab=pc1lab, ylab=pc2lab, main=main, ...)

	if (plotlabels){
		# Add text to the points
		with(as.data.frame(pca$x), textxy(PC1, PC2, labs=nnum, cex=0.7))

		# setup for no margins on the legend
		# c(bottom, left, top, right)
		par(mar=c(0, 0, 0, 0))
		plot.new()

		legend('center','groups', paste(nnum, rnames, sep=" = "), xpd=TRUE, ncol=2, cex=0.5, text.col=colors[fac], bty = "n")
		legend("top", legend=levels(fac), col=colors, pch=20, xpd=TRUE, ncol=2, cex=0.5, bty = "n")

		# Restore default clipping rect
		par(mar=c(5, 4, 4, 2) + 0.1)
	} else{
		# Expand right side of clipping rect to make room for the legend
		legend("top",                   # Location of legend 
			xpd = TRUE,                      # Allow drawing outside plot area
			xjust = 0,                       # Left justify legend box on x
			yjust = .5,                      # Center legend box on y
			legend = levels(fac),            # legend element labels
			col = colors,                    # Legend Element colors
			pch = 20,                        # Legend Element Styles
			cex = 0.5,
			ncol=2
		)

	}

	par(new=TRUE)
	# Turn off device driver (to flush output to PNG file)
	dev.off()
}

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="top", labelsig=TRUE, textcx=1, ...) {
	with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
	with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
	with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
	with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
	if (labelsig) {
		require(calibrate)
		with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
	}
	legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR < ",sigthresh,sep=""), paste("abs(LogFC) > ",lfcthresh,sep=""), "both "), pch=20, col=c("red","orange","green"), ncol=3)

	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_volcanoPlot <- function(resdata, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Volcano Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")
	
	# Plot scaled volcano plots 
	# Create the png filename and start PNG device driver to save output to figure.png
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_volcanoplot.png",sep='');
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), ylim=c(0, 15), labelsig=FALSE, main=mtitle)
	
	# Create the pdf filename and start PDF device driver to save output to figure.pdf
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_volcanoplot.pdf",sep='')
	pdf(pdffile)
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), ylim=c(0, 15), labelsig=FALSE, main=mtitle)

	# Plot auto scaled volcano plots 
	# Create the png filename and start PNG device driver to save output to figure.png
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_autoscaled_volcanoplot.png",sep='');
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, labelsig=FALSE, main=mtitle)
	
	# Create the pdf filename and start PDF device driver to save output to figure.pdf
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_autoscaled_volcanoplot.pdf",sep='')
	pdf(pdffile)
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, labelsig=FALSE, main=mtitle)
}

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
maplot <- function (res, main="MA Plot", thresh=0.05, labelsig=TRUE, textcx=1, ...) {
	with(res, plot(baseMean + 0.1, log2FoldChange, pch=20, cex=.5, log="x", main=main, ...))
	with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
	if (labelsig) {
		require(calibrate)
		with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
	}
	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_maPlot <- function(resdata, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("MA Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")

	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_maplot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	maplot(resdata, main=mtitle, labelsig=FALSE)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_maplot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	maplot(resdata, main=mtitle, labelsig=FALSE)
}

## Examine plot of p-values
pvalplot <- function (res, main="p-values", ...) {
	#hist(res$pvalue, breaks=50, col="grey", main=main)
	hist(res$pvalue[res$baseMean > 1], breaks=50, col="grey", main=main)
	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_pvalPlot <- function(res, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Pvalues Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")
	
	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_pvalplot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	pvalplot(res, main=mtitle)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_pvalplot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	pvalplot(res, main=mtitle)
}

## Examine independent filtering
indfilplot <- function (res, main="Independent filtering",...) {
	plot(metadata(res)$filterNumRej, type="b", xlab="quantiles of 'baseMean'", ylab="number of rejections", main=main)
	par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_indfilPlot <- function(res, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Independent Filtering Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")

	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_independent_filtering_plot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=8, res=300, units="in");
	# Plot PCA without labels
	indfilplot(res, main=mtitle)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_independent_filtering_plot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	indfilplot(res, main=mtitle)
}

## Heatmap of sample-to-sample distances using the Poisson Distance
sampleDistplot <- function (fdds, rld, main="Independent filtering",...) {
	# Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of 
	# dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples.
	# The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns, so we need to transpose the 
	# counts in dds.
	suppressMessages(library("PoiClaClu"))
	suppressMessages(library("pheatmap"))
	poisd <- PoissonDistance(t(counts(fdds)))

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
         fontsize=8,
         main=main)

	#par(new=TRUE)
	# Turn off device driver (to flush output to PNG/PDF file)
	dev.off()
}

plot_sampleDistplot <- function(fdds, rld, digPlotDirPng, digPlotDirPdf, outputFile){
	mtitle <- paste("Sample Correlation (Poisson Distance) Plot (", conditionVals[2],", ", conditionVals[1],")", sep="")

	# Create the png filename
	pngfile  <- paste(digPlotDirPng, "/", basename(outputFile), "_sampleCorrelation_poissonDistance_plot.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(filename=pngfile, height=8,width=16, res=300, units="in");
	# Plot PCA without labels
	sampleDistplot(fdds, rld, main=mtitle)

	# Create the pdf filename
	pdffile<-paste(digPlotDirPdf, "/", basename(outputFile), "_sampleCorrelation_poissonDistance_plot.pdf",sep='')
	# Start PDF device driver to save output to figure.pdf
	pdf(pdffile)
	# Plot PCA with labels
	sampleDistplot(fdds, rld, main=mtitle)
}


