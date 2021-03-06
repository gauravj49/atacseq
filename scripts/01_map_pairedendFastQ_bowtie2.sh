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

# Expected naming convention: ExperimentInfo_Species_SeqProtocol_SeqType_Read
# 1) Experiment info: ExperimentName-CellLine-SpecialCondition-ReplicateNumber
# 2) Species: Human, mouse, etc.
# 3) SeqProtocol: Rnaseq, ChIPseq etc.
# 4) SeqType: paired end or single end
# 5) Read: R1 for side 1 and R2 for side 2
# Example: 
#   - TransPB-CD4-CD4PosTcells-Rep3_mm_atacseq_se_R1.fastq.gz
#   - TransPB-CD4-GSM2056292-Rep1_mm_atacseq_pe_R1.fastq.gz and TransPB-CD4-GSM2056292-Rep1_mm_atacseq_pe_R2.fastq.gz

for f in ${fastqdir}/*_R1.fastq.gz
do 
 # Get basename name
 bname=$(basename "${f}" .fastq.gz)
 fq1=${f}
 fq2=${fq1/_R1/_R2}
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

 echo "# Trim the adapters using TrimGalore" >> "${scriptFile}"
 echo "trim_galore -j 4 --fastqc --fastqc_args \"--outDir ${trimFastqcDir} \" -o ${trimmedFastDir} --paired ${fq1} ${fq2} 2>&1 | tee ${mappingLogFile}" >> "${scriptFile}"
 echo "" >> "${scriptFile}"

 # Rename the fastq files that has *_val_1.fq.gz and *_val_2.fq.gz
 echo "# Rename the fastq files that has *_val_1.fq.gz and *_val_2.fq.gz" >> "${scriptFile}"
 echo "mv ${trimmedFastDir}/${bname}_val_1.fq.gz ${trimmedFastDir}/${bname}.fastq.gz 2>&1 | tee -a ${mappingLogFile}" >> "${scriptFile}" 
 echo "mv ${trimmedFastDir}/$(basename "${f}" _1.fastq.gz)_2_val_2.fq.gz ${trimmedFastDir}/$(basename "${f}" _1.fastq.gz)_2.fastq.gz 2>&1 | tee -a ${mappingLogFile}" >> "${scriptFile}" 
 echo "" >> "${scriptFile}"

 # Align reads
	 #$BT2_HOME/bowtie2 --local -x specious -U $BT2_HOME/example/reads/longreads.fq -S eg.sam
	 #-U Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. Reads may be a mix of different lengths. If - is specified, bowtie2 gets the reads from the "standard in" or "stdin" filehandle.
	 #-x The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.
	 #-S File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).
 echo "# Align reads" >> "${scriptFile}"
 echo "bowtie2 -x ${bt2index} -p 8 --very-sensitive -1 ${fq1} -2 ${fq2} -S ${bamDir}/${bname}.sam  2>&1 | tee -a ${mappingLogFile}" >> "${scriptFile}"
 echo "" >> "${scriptFile}"

 # Convert the SAM to BAM
 echo "# Convert the SAM to BAM" >> "${scriptFile}"
 echo "samtools view -b -S ${bamDir}/${bname}.sam > ${bamfile}.unsorted" >> "${scriptFile}"
 echo "" >> "${scriptFile}"
 
 # Sort the bam file
 echo "# Sort the bam file" >> "${scriptFile}"
 echo "samtools sort -o ${bamfile} ${bamfile}.unsorted" >> "${scriptFile}"
 echo "" >> "${scriptFile}"
 
 # Remove extra files
 echo "# Remove extra files" >> "${scriptFile}"
 echo "rm -rf ${bamDir}/${bname}.sam ${bamfile}.unsorted"  >> "${scriptFile}"
 echo "" >> "${scriptFile}"

 # Remove duplicates from the bam file
 echo "# Remove duplicates from the bam file" >> "${scriptFile}"
 echo "samtools rmdup -s ${bamfile} ${rmdupbamfile}.unsorted" >> "${scriptFile}"
 echo "" >> "${scriptFile}"

 # Sort the bam file
 echo "# Sort the bam file" >> "${scriptFile}"
 echo "samtools sort -o ${rmdupbamfile} ${rmdupbamfile}.unsorted" >> "${scriptFile}"
 echo "" >> "${scriptFile}"

 # Remove extra files again
 echo "# Remove extra files again" >> "${scriptFile}"
 echo "rm -rf ${bamfile} ${rmdupbamfile}.unsorted"  >> "${scriptFile}"
 echo "" >> "${scriptFile}"

 # create index of the bam files
 echo "# create index of the bam files" >> "${scriptFile}"
 echo "samtools index ${rmdupbamfile}" >> "${scriptFile}"
 echo "" >> "${scriptFile}"

 # Write the command in the script file and give it correct permission to run
 chmod 775 "${scriptFile}"
 echo "${scriptFile}"

done

