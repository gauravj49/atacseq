#! /usr/bin/env Rscript

# ****************************************************************
# GOAL   : Identify COREs using CREAM (Clustering of genomic REgions Analysis Method)
# USAGE  : Rscript scripts/R_identify_COREs_using_CREAM.R -if=cream/individual/peaks/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup_peaks.narrowPeak -of=cream/individual/cores/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup_peaks_COREs.bed
# SOURCE : https://github.com/bhklab/CREAM
# ****************************************************************

################ load libraries ###############
suppressPackageStartupMessages(library("argparse"))

## Main logic
main <- function(){    

  ##### Parse command line arguments ############
  args        <- check_options()
  inputfile   <- args$inputfile
  outputfile  <- args$outputfile
  MinLength   <- args$MinLength
  peakNumMin  <- args$peakNumMin
  WScutoff    <- args$WScutoff

  cat("\n- Input arguments")
  cat("\n\t- inputfile  = ", inputfile)
  cat("\n\t- outputfile = ", outputfile)
  cat("\n\t- MinLength  = ", MinLength)
  cat("\n\t- peakNumMin = ", peakNumMin)
  cat("\n\t- WScutoff   = ", WScutoff)
  ################ Load libraries ###############
  load_libraries()

  ################ Main logic ###################
  ## Get the input data
  cat("\n\n- Reading input file ...\n")
  options(warn=-1) # Turn off warnings
  IdentifiedCOREs <- CREAM(in_path = inputfile, MinLength = 1000, peakNumMin = 2, WScutoff = 0.25)
  options(warn=0)  # Turn on  warnings

  # Save the data frame into a output file
  coreDF           <- data.frame(IdentifiedCOREs[,c(1:3,6)])
  colnames(coreDF) <- c('chr','start','end','COREs_score')
  
  # Add a name column
  coreDF$name      <- paste(coreDF$chr, coreDF$start, coreDF$end, coreDF$COREs_score, sep='_')

  # Reorder the columns
  coreDF           <- coreDF[, c('chr','start','end','name', 'COREs_score')]

  # Write to the output file
  fwrite(coreDF, file=outputfile, quote=F, row.names=F, col.names=F, sep='\t',  nThread=48)

  # Sort file
  system(paste0("sort -k1,1 -k2n ", outputfile, " -o ", outputfile))

  # # Print session info to the session's log file
  # logdir <- paste(dirname(outputfile),"/logs", sep='')
  # system(paste("mkdir -p ", logdir, sep=''))
  # session_logfile <- paste(logdir,"/session_info_",basename(outputfile),".log" ,sep='');
  # print_session_info(session_logfile)
}

##################### USER DEFINIED FUNCTIONS #########################

load_libraries <- function(){
	# Load libraries at start
suppressPackageStartupMessages(library(CREAM))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(data.table))

}

# 
# https://community.rstudio.com/t/best-practices-for-saving-session-information/44836/2

# Print session info as log file formatted in tabular format
# Source: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
print_session_info_old <- function(session_logfile){
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
		- Rscript scripts/R_identify_COREs_using_CREAM.R -if=cream/individual/peaks/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup_peaks.narrowPeak -of=cream/individual/cores/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup_peaks_COREs.bed
	-------------------------------------------------
	CONTACT: 
		Gaurav Jain
		gaurav.jain@tum.de
	-------------------------------------------------\n
	")
	# create parser object
	parser <- ArgumentParser(description=cat(desc))

	# Add arguments 
	parser$add_argument("-if", "--inputfile"  , dest="inputfile"  , help="*Input counts matrix file", type="character", required=TRUE)
	parser$add_argument("-of", "--outputfile" , dest="outputfile" , help="*Output file name" , type="character", required=TRUE)
  parser$add_argument('-ml', "--minLength"  , dest="MinLength"  , help=" Criteria for the minimum number of functional regions in the input file (Default=%(default)s)", default=1000)
  parser$add_argument('-pm', "--peakNumMin" , dest="peakNumMin" , help=" Minimum number of peaks for CORE identification (Default=%(default)s)", default=2)
  parser$add_argument('-wc', "--WScutoff"   , dest="WScutoff"   , help=" Threshold used to identify WS within distribution of maximum distance between peaks for each order of CORE (Default=%(default)s)", default=0.25)
	

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

