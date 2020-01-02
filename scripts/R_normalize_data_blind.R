#! /usr/bin/env Rscript

# ****************************************************************
# GOAL  : Normalize with various methods:
# 		  - Quantile normalization
# USAGE : 
# ****************************************************************

################ load libraries ###############
suppressPackageStartupMessages(library("argparse"))

## Main logic
main <- function(){    

	##### Parse command line arguments ############
	args        <- check_options()
	inputfile   <- args$inputfile
	outputfile  <- args$outputfile

	################ Load libraries ###############
	load_libraries()

	################ Main logic ###################
	## Get the input data
	cat("- Reading input file ...\n")
	allDataOrig     <- data.frame(read.table(inputfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))

	# Create the output directories and files
	cat("\n2) Create the output directories and files ...\n")
	outputDir       <- paste(dirname(outputFile),"/snormalized_counts", sep=''); system(paste("mkdir -p", outputDir, sep=' '))
	basefilename    <- tools::file_path_sans_ext(basename(outputFile))
	extfilename     <- tools::file_ext(outputFile)

	# Print session info to the session's log file
	logdir <- paste(dirname(outputfile),"/logs", sep='')
	system(paste("mkdir -p ", logdir, sep=''))
	session_logfile <- paste(logdir,"/session_info_",file_path_sans_ext(basename(outputfile)),".log" ,sep='');
	print_session_info(session_logfile)
}

############ USER DEFINED FUNCTIONS ##########
# Generate RLE plots 
get_rle_plots <- function(sampleCountsDF, outputDir, dataType, color){
	# # Relative log expression (RLE) plot from ExpressionSet:
	# 	- RLE plots were initially proposed to measure the overall quality of a dataset but can also be used to visualize the presence of unwanted batch effects in the data
	#  	- Unwanted variation can be highly problematic and so its detection is often crucial
	#  	- Relative log expression (RLE) plots are a powerful tool for visualising such variation in high dimensional data
	#  	- RLE plots are particularly useful for assessing whether a procedure aimed at removing unwanted variation, i.e. a normalisation procedure, has been successful
	#  	- These plots, while originally devised for gene expression data from microarrays, can also be used to reveal unwanted variation in single-cell expression data, where such variation can be problematic

	# # If style is "full", as usual with boxplots:
	# 	- The box shows the inter-quartile range and whiskers extend no more than 1.5 * IQR from the hinge (the 25th or 75th percentile)
	# 	- Data beyond the whiskers are called outliers and are plotted individually
	# 	- The median (50th percentile) is shown with a white bar

	# Get the output dirs
	rlePlotsDir     <- paste(outputDir, "/rle_plots", sep=''); system(paste("mkdir -p", rlePlotsDir, sep=' '));
	ylimvals        <- c(-5,5)

	# 1) # Create RLE plots for RAW counts
	cat("\t9.1) Create RLE plots for RAW counts ...\n")
	filter          <- apply(sampleCountsDF, 1, function(x) length(x[x>1])>=2)
	filtered        <- sampleCountsDF[filter,]
	set             <- newSeqExpressionSet(as.matrix(filtered))
	RLEplotFile     <- paste(rlePlotsDir, "/", basefilename, "_",dataType,".png", sep='') 
	png(filename=RLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(set, outline=FALSE,col=color, las=2, cex.axis=0.5, ylim=ylimvals)
	dev.off()

}


# Get normalized reads 
get_qn_reads <- function(fdds, sampleCountsDF, outputFile){
	# 1) RLE plots for raw data
	# 2) Get size factor normalized reads (SFN)
	# 3) Get variance stabilizing transformation (VST) normalized reads

	# Get annotation table
	sampleFiles     <- c(ControlFiles, TreatmentFiles)
	sampleCondition <- c(rep(conditionVals[1], length(ControlFiles)), rep(conditionVals[2], length(TreatmentFiles)))
	sampleDiagnosis <- c(rep(0, length(ControlFiles)), rep(1, length(TreatmentFiles)))
	newSampleTable  <- data.frame(sampleName = names(sampleCountsDF), condition = sampleCondition, diagnosis = sampleDiagnosis, row.names=names(sampleCountsDF))

	# 1) # Create RLE plots for RAW counts
	cat("\t9.1) Create RLE plots for RAW counts ...\n")
	filter          <- apply(sampleCountsDF, 1, function(x) length(x[x>1])>=2)
	filtered        <- sampleCountsDF[filter,]
	set             <- newSeqExpressionSet(as.matrix(filtered))
	rawRLEplotFile  <- paste(rlePlotsDir, "/", basefilename, "_raw.png", sep='') 
	png(filename=rawRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(set, outline=FALSE,col='navy',las=2, cex.axis=0.5, ylim=ylimvals)
	dev.off()

	# 2) Get quantile normalized reads (QNR)
	cat("\t9.2) Get quantile normalized reads (QNR) ...\n")
	suppressPackageStartupMessages(library("preprocessCore", warn.conflicts=FALSE, quietly=TRUE))
	qnrCountsDF       <- normalize.quantiles(as.matrix(sampleCountsDF))
	qnrfilter         <- apply(qnrCountsDF, 1, function(x) length(x[x>1])>=2)
	qnrfiltered       <- qnrCountsDF[qnrfilter,]
	qnrset            <- newSeqExpressionSet(as.matrix(qnrfiltered))

	qnrrawRLEplotFile <- paste(rlePlotsDir, "/", basefilename, "_qnr.png", sep='') 
	png(filename=qnrrawRLEplotFile, height=900, width=1000, bg="white", res=100)
	par(mar=c(20,4,1,1));
	plotRLE(qnrset, outline=FALSE,col='deepskyblue',las=2, cex.axis=0.5, ylim=ylimvals)
	dev.off()

	# Get QNR related directories and filenames
	qnrnormalizedDir      <- paste(outputDir,'/qnr/',conditionVals[2],'_',conditionVals[1], sep='')
	qnrtnormalizedDir     <- paste(qnrnormalizedDir, "/transposed", sep='')
	indqnrnormalizedDir   <- paste(qnrnormalizedDir, "/qnr_normalized_counts", sep='')
	qnrnormalizedOutFile  <- paste(qnrnormalizedDir , "/", basefilename, "_qnr_normalized.txt", sep='')
	qnrtnormalizedOutFile <- paste(qnrtnormalizedDir, "/", basefilename, "_featureInCols_qnr_normalized.txt", sep='')
	system(paste("mkdir -p", qnrnormalizedDir, qnrtnormalizedDir, indqnrnormalizedDir, sep=' '))

	# Convert rownames as feature column to save it in the csv file
	qnrrownames <- rownames(sampleCountsDF[qnrfilter,])
	qnrcolnames <- names(sampleCountsDF[qnrfilter,])
	colnames(qnrfiltered) <- qnrcolnames
	rownames(qnrfiltered) <- qnrrownames
	qnrfiltered <- cbind(feature = rownames(qnrfiltered), qnrfiltered)
	rownames(qnrfiltered) <- 1:nrow(qnrfiltered)
	
	# Save the counts file to output csv file
	write.table(qnrfiltered   , file = qnrnormalizedOutFile , row.names = F, sep = '\t', quote = F)
	write.table(t(qnrfiltered), file = qnrtnormalizedOutFile, col.names = F, sep = '\t', quote = F)

	# Get the normalized individual count files
	indvcmd <- paste("bash scripts/get_individual_countFiles.sh", qnrnormalizedOutFile, indqnrnormalizedDir, rnaClass, 1, 2, sep=' ')
	system(indvcmd)

}

get_counts_matrix <- function(inputDir, rnaClass){
	# Combine individual HTSeq count files into a matrix
	# Source: https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4

	# Get the list of all count files
	allCountFiles <- list.files(inputDir)

	# Subset the files here for specific files
	subCountFiles <- allCountFiles # replace it with the subset file code

	# Read each file as array element of DT and rename the last 2 cols 
	# we created a list of single sample tables
	DT <- list()
	for (i in 1:length(subCountFiles)){
		# Get the absolute file name
		infile = paste(inputDir, subCountFiles[i], sep = "/")

		# Add it to the data table
		DT[[subCountFiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)

		# Get the column names 
		colnames(DT[[subCountFiles[i]]]) <- c("feature", subCountFiles[i])
	}

	# Merge all elements based on first 'feature' columns
	mergedDataTab <- DT[[subCountFiles[1]]]
 
	# We now add each other table with the 'feature' column as key
	for (i in 2:length(subCountFiles)) {
		y <- DT[[subCountFiles[i]]]
		z <- merge(mergedDataTab, y, by = c("feature"))
		mergedDataTab <- z
	}
 
	# 'feature' column becomes rownames
	rownames(mergedDataTab) <- mergedDataTab$feature
	mergedDataTab <- mergedDataTab[,-1]
	 
	## add total counts per sample
	mergedDataTab <- rbind(mergedDataTab, tot.counts=colSums(mergedDataTab))

	# Remove rnaClass (mirna or allncrna) from column names
	oldColNames            <- names(mergedDataTab)
	newColNames			   <- gsub(paste("_",rnaClass,"Counts.txt", sep=''), "" , oldColNames)
	names(mergedDataTab)  <- newColNames

	# Remove HTseq specific rows
	# __no_feature,__not_aligned,__too_low_Qual,__alignment_not_unique,__ambiguous
	toDrop <- c('tot.counts','__no_feature','__not_aligned','__too_low_Qual','__alignment_not_unique','__ambiguous')
	mergedDataTab <- mergedDataTab[ !(rownames(mergedDataTab) %in% toDrop), ]

	return(mergedDataTab)

}
##################### USER DEFINIED FUNCTIONS #########################
# Load Libraries
load_libraries <- function(){
	# Load libraries at start
	suppressPackageStartupMessages(library("preprocessCore", warn.conflicts=FALSE, quietly=TRUE))
	suppressPackageStartupMessages(library("edgeR" , warn.conflicts=FALSE, quietly=TRUE))
	suppressPackageStartupMessages(library("RUVSeq", warn.conflicts=FALSE, quietly=TRUE))
	# suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))

}

# Print session info as log file formatted in tabular format
# Source: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
print_session_info <- function(session_logfile){
	suppressPackageStartupMessages(library("devtools"))
	suppressPackageStartupMessages(library("knitr"))

	# Get all the session info to the variable
	my_session_info <- devtools::session_info()

	# Print it in the tabular format using knitr
	writeLines(text = {
	    paste(sep = "\n", collapse = "",
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 32), collapse = ""),
	                "R environment",
	                paste0(rep("-", 33), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(data.frame(setting = names(my_session_info$platform),
	                                  value = as.character(my_session_info$platform))), collapse = "\n"),
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 35), collapse = ""),
	                "packages",
	                paste0(rep("-", 35), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(my_session_info$packages), collapse = "\n")
	          )
	}, con = session_logfile)
}

check_options <- function(){
	# Description of the script
	desc <- sprintf("
	----------------- SAMPLE USAGE ------------------
		- Rscript scripts/R_heatmap.R -if=output/filtered_data/deseq2/mirna/results_DEseq2/ACC18mwt_over_ACC4mwt_DE_RESULTS.txt -of=output/filtered_data/deseq2/mirna/heatmaps_DEseq2/ACC18mwt_over_ACC4mwt_DE_heatmaps.txt -xc -yc
		- Rscript scripts/R_heatmap.R -if=output/filtered_data/deseq2/mirna/results_DEseq2/Blood18mwt_over_Blood4mwt_DE_RESULTS.txt -of=output/filtered_data/deseq2/mirna/heatmaps_DEseq2/Blood18mwt_over_Blood4mwt_DE_heatmaps.txt
	-------------------------------------------------
	CONTACT: 
		Gaurav Jain
		gaurav.jain@dzne.de
	-------------------------------------------------\n
	")
	# create parser object
	parser <- ArgumentParser(description=cat(desc))

	# Add arguments 
	parser$add_argument("-if", "--inputfile"  , dest="inputfile"  , help="*Input counts matrix file", type="character", required=TRUE)
	parser$add_argument("-of", "--outputfile" , dest="outputfile" , help="*Output file name" , type="character", required=TRUE)
	

	# Print the help message only if no arguments are supplied
	if(length(commandArgs(TRUE))==0){
		cat(desc)
		parser$print_help()
		quit()
	}

	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	args <- parser$parse_args()
	return(args)
}

## Call the main function in the end
main()

