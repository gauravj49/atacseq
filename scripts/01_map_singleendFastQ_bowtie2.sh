#!/bin/bash

# USAGE: bash 01_map_pairedendFastQ_bowtie2.sh input/fastq/hgStomachF35 output/hgStomachF35 mm10

# Set user defined environment variables
jobdir="/home/rad/users/gaurav/projects/seqAnalysis/atacseq"
fastqdir=$1
outDir=$2
projName=$3
species=$4
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/${projName}"
bt2index="/home/rad/packages/bowtie2/indexes/${species}"
bamDir="${outDir}/bams/trimmed"
origFastqcDir="${outDir}/qc/fastqc/00_original"
trimFastqcDir="${outDir}/qc/fastqc/01_trimmed"
mappingLogsDir="${outDir}/logs/trimmingLogs"
trimmedFastDir="${outDir}/tempLocal/trimmed_fastq"

# echo "${scriptsdir} ${bamDir} ${mappingLogsDir} ${trimmedFastDir} ${origFastqcDir} ${trimFastqcDir}"

# Create required dirs
mkdir -p ${scriptsdir} ${bamDir} ${mappingLogsDir} ${trimmedFastDir} ${origFastqcDir} ${trimFastqcDir}

for f in ${fastqdir}/*.fastq.gz
do 
 # Get basename name
 bname=$(basename "$f" .fastq.gz)
 fq1=${f}
 scriptFile="${scriptsdir}/${bname}.sh"
 bamfile="${bamDir}/${bname}.bam"
 rmdupbamfile="${bamDir}/${bname}_rmdup.bam"
 mappingLogFile="${mappingLogsDir}/${bname}_trimming.log"

 # Get the jobname to submit for each job
 jobname="01_$bname"

 # Create the script file
 touch "${scriptFile}"
 echo "#!/bin/bash" > "${scriptFile}"
 echo "" >> "${scriptFile}"

 # Trim the adapters using TrimGalore
 #  -j/--cores INT:  It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.
 #  -o/--output_dir <DIR>: If specified all output will be written to this directory instead of the current directory. If the directory doesn't exist it will be created for you.
 # --fastqc : Run FastQC in the default mode on the FastQ file once trimming is complete.
 # --fastqc_args "<ARGS>": Passes extra arguments to FastQC. If more than one argument is to be passed to FastQC they must be in the form arg1 arg2 [..].
 #   					   An example would be: --fastqc_args "--nogroup --outDir /home/".

 echo "trim_galore -j 4 --fastqc --fastqc_args \"--outDir ${trimFastqcDir} \" -o ${trimmedFastDir} ${f} 2>&1 | tee ${mappingLogFile}" >> "${scriptFile}"

 # Align reads
	 #$BT2_HOME/bowtie2 --local -x specious -U $BT2_HOME/example/reads/longreads.fq -S eg.sam
	 #-U Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. Reads may be a mix of different lengths. If - is specified, bowtie2 gets the reads from the "standard in" or "stdin" filehandle.
	 #-x The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.
	 #-S File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).
 echo "bowtie2 -x ${bt2index} -p 8 --very-sensitive -U ${fq1} > ${bamDir}/${bname}.sam  2>&1 | tee -a ${mappingLogFile}" >> "${scriptFile}"

 # Convert the SAM to BAM
 echo "samtools view -b -S ${bamDir}/${bname}.sam > ${bamfile}.unsorted" >> "${scriptFile}"

 # Sort the bam file
 echo "samtools sort -o ${bamfile} ${bamfile}.unsorted" >> "${scriptFile}"

 # Remove extra files
 echo "rm -rf ${bamDir}/${bname}.sam ${bamfile}.unsorted"  >> "${scriptFile}"

 # Remove duplicates from the bam file
 echo "samtools rmdup -s ${bamfile} ${rmdupbamfile}.unsorted" >> "${scriptFile}"

 # Sort the bam file
 echo "samtools sort -o ${rmdupbamfile} ${rmdupbamfile}.unsorted" >> "${scriptFile}"

 # Remove extra files again
 echo "rm -rf ${bamfile} ${rmdupbamfile}.unsorted"  >> "${scriptFile}"

 # create index of the bam files
 echo "samtools index ${rmdupbamfile}" >> "${scriptFile}"

 # Write the command in the script file and give it correct permission to run
 chmod 775 "${scriptFile}"
 echo "${scriptFile}"

done

