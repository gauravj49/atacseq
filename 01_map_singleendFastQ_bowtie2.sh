#!/bin/bash

# USAGE: bash 01_map_pairedendFastQ_bowtie2.sh input/fastq/hgStomachF35 output/hgStomachF35 hg19

# Set user defined environment variables
jobdir="/usr/users/gjain/bin/projects/atacseq"
fastqdir=${jobdir}/$1
bamdir=${jobdir}/$2
species=$3
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2"
stdoutFile="${bamdir}/map_fastQ_bowtie2.out" 
errorsFile="${bamdir}/map_fastQ_bowtie2.err" 
bt2index="/usr/users/gjain/bin/libs/bowtie2/indexes/${species}/${species}"

# Create required dirs
mkdir -p ${scriptsdir} ${bamdir}

for f in ${fastqdir}/*_R1_*.fastq.gz
do 
 # Get bamfile name
 fastQfile=$(basename "$f" .fastq.gz)
 fq1=${f}
 scriptFile="${scriptsdir}/${fastQfile}.sh"
 bamfile=${bamdir}/${fastQfile}.bam

 # Get the jobname to submit for each job
 jobname="01_$fastQfile"

 # Create the script file
 touch "${scriptFile}"

 # Align reads
	 #$BT2_HOME/bowtie2 --local -x specious -U $BT2_HOME/example/reads/longreads.fq -S eg.sam
	 #-U Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. Reads may be a mix of different lengths. If - is specified, bowtie2 gets the reads from the "standard in" or "stdin" filehandle.
	 #-x The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.
	 #-S File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).
 echo "bowtie2 -x ${bt2index} -p 8 --very-sensitive -U ${fq1} > ${bamdir}/${fastQfile}.sam" > "${scriptFile}"

 # Convert the SAM to BAM
 echo "samtools view -b -S ${bamdir}/${fastQfile}.sam > ${bamdir}/${fastQfile}.unsorted.bam" >> "${scriptFile}"

 # Sort the bam file
 echo "samtools sort -o ${bamdir}/${fastQfile}.bam ${bamdir}/${fastQfile}.unsorted.bam" >> "${scriptFile}"

 # Remove extra files
 echo "rm -rf ${bamdir}/${fastQfile}.sam ${bamdir}/${fastQfile}.unsorted.bam"  >> "${scriptFile}"

 # create index of the bam files
 echo "samtools index ${bamfile} " >> "${scriptFile}"

 # Write the command in the script file and give it correct permission to run
 chmod 775 "${scriptFile}"

 # Submit the job
 echo -e "Submitting job for ${bamfile}" 
 bsub -q fat -W 36:00 -n 8 -e "${errorsFile}" -o "${stdoutFile}" -J ${jobname} sh "${scriptFile}"
 echo -e " "
done
