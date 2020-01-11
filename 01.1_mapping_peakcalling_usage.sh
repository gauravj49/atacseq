# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

############### For Users ###################
# Get the wrappers and then run them in parallel

# 1) Christine Klement
# 1.1) For the first run for ATACseq_MUC001 sampels
species="mm10"
user="christine"
projName="AGRad_ATACseq_MUC001"
fastqdir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/demultiplexing/output/fastq"
outputDir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001"
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/${projName}"
multiqcDir="${outputDir}/qc/multiqc"; mkdir -p ${multiqcDir}
fastqcDir="${outputDir}/qc/fastqc"; mkdir -p ${fastqcDir}
logsDir="${outputDir}/logs"; mkdir -p ${logsDir}

bash scripts/01_map_singleendFastQ_bowtie2.sh ${fastqdir} ${outputDir} ${projName} ${species}
cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Run multiqc on original files
ls ${fastqdir}/*.gz | parallel --progress --eta -j 64 'fastqc -o ${fastqcDir}/00_original {}'
multiqc -o ${multiqcDir} -n 01_original_fastq_stats_${projName} ${fastqcDir}/00_original

# Run multiqc on trimmed files
multiqc -o ${multiqcDir} -n 02_trimmed_fastq_stats_${projName} ${fastqcDir}/01_trimmed ${logsDir}/trimmingLogs

# Get mapping stats and run multiqc on mapping logs
bamdir="${outputDir}/bams/trimmed"
ls ${bamdir}/*.bam | parallel --progress --eta -j 16 'samtools stats -@ 32 {} > {.}_mapping.logs'
mappingLogsDir="${outputDir}/logs/mappingLogs"; mkdir -p ${mappingLogsDir}; mv ${bamdir}/*.logs ${mappingLogsDir}
multiqc -o ${multiqcDir} -n 03_mapping_stats_${projName} ${mappingLogsDir}

# Get peaks and homer annotation
jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="mm10"
user="christine"
projName="AGRad_ATACseq_MUC001"
fastqdir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/demultiplexing/output/fastq"
outputDir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001"
bamdir="${outputDir}/bams/trimmed"
peaksdir="${outputDir}/peaks"
scriptsdir="${jobdir}/scripts/02_peakCallingMACs_annotationHomer/${projName}"
mkdir -p ${scriptsdir}
for b in ${bamdir}/*.bam; do bash scripts/02_peakCallingMACS2_annotationHomer.sh ${b} ${peaksdir} ${species} ${projName}; done;
cmd="parallel --tmpdir /media/rad/SSD1/atac_temp ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}


# 2) Anja Pfaus
# 2.1) For T-ALL human cell lines
species="hg19"
user="anja"
projName="tALLcellLine_hg"
fastqdir="/media/rad/HDD1/PUB_CRCs/atac/atac_tALLcellLine_hg"
outputDir="/media/rad/HDD1/atacseq/output/anja/tALLcellLine_hg"
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/${projName}"
multiqcDir="${outputDir}/qc/multiqc"; mkdir -p ${multiqcDir}

bash scripts/01_map_singleendFastQ_bowtie2.sh ${fastqdir} ${outputDir} ${projName} ${species}
cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Run multiqc on trimmed files
trimFastqcDir="${outputDir}/qc/fastqc/01_trimmed"
multiqc -o ${multiqcDir} -n 01_trimmed_fastq_stats_${projName} ${trimFastqcDir}

# Run multiqc on trimming logs 
trimmingLogsDir="${outputDir}/logs/trimmingLogs"
multiqc -o ${multiqcDir} -n 02_trimming_stats_${projName} ${trimmingLogsDir}

# Get mapping stats and run multiqc on mapping logs
bamdir="${outputDir}/bams/trimmed"
ls ${bamdir}/*.bam | parallel --progress --eta -j 16 'samtools stats -@ 32 {} > {.}_mapping.logs'
mappingLogsDir="${outputDir}/logs/mappingLogs"; mkdir -p ${mappingLogsDir}; mv ${bamdir}/*.logs ${mappingLogsDir}
multiqc -o ${multiqcDir} -n 03_mapping_stats_${projName} ${mappingLogsDir}

# Get peaks and homer annotation
jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="hg19"
projName="tALLcellLine_hg"
fastqdir="/media/rad/HDD1/PUB_CRCs/atac/atac_tALLcellLine_hg"
outputDir="/media/rad/HDD1/atacseq/output/anja/tALLcellLine_hg"
bamdir="${outputDir}/bams/trimmed"
peaksdir="${outputDir}/peaks"
scriptsdir="${jobdir}/scripts/02_peakCallingMACs_annotationHomer/${projName}"
mkdir -p ${scriptsdir}
for b in ${bamdir}/*.bam; do bash scripts/02_peakCallingMACS2_annotationHomer.sh ${b} ${peaksdir} ${species} ${projName}; done;
cmd="parallel --tmpdir /media/rad/SSD1/atac_temp ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# 2.2) For T-ALL mouse cell lines
jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="mm10"
user="sabrina"
projName="nkTimecourse"
projDir="/home/rad/media/rad/HDD1/atacseq/sabrina/nkTimecourse"
fastqdir="${projDir}/fastq"
mappingDir="${projDir}/mapping"
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/${projName}"
multiqcDir="${projDir}/qc/multiqc"; 

# Create relevant dirs
mkdir -p ${mappingDir} ${scriptsdir} ${multiqcDir}

# Perform mapping
bash scripts/01_map_singleendFastQ_bowtie2.sh ${fastqdir} ${mappingDir} ${projName} ${species}
cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}


# 3) Sabrin Bortoluzzi from AG Schmidt-Supprian 
# 3.1) For the NK-Timecouse ATACseq samples
jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="mm10"
user="sabrina"
projName="nkTimecourse"
projDir="/home/rad/media/rad/HDD1/atacseq/sabrina/nkTimecourse"
fastqdir="${projDir}/fastq"
mappingDir="${projDir}/mapping"
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/${projName}"
multiqcDir="${projDir}/qc/multiqc"; 

# Create relevant dirs
mkdir -p ${mappingDir} ${scriptsdir} ${multiqcDir}

# Perform mapping
bash scripts/01_map_singleendFastQ_bowtie2.sh ${fastqdir} ${mappingDir} ${projName} ${species}
cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Get peaks and homer annotation
jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="mm10"
user="sabrina"
projName="nkTimecourse"
projDir="/home/rad/media/rad/HDD1/atacseq/sabrina/nkTimecourse"
fastqdir="${projDir}/fastq"
mappingDir="${projDir}/mapping"
bamdir="${projDir}/mapping/bams/trimmed"
peaksdir="${projDir}/peaks"
scriptsdir="${jobdir}/scripts/02_peakCallingMACs_annotationHomer/${projName}"
mkdir -p ${scriptsdir} ${bam} ${peaksdir}
for b in ${bamdir}/*.bam; do bash scripts/02_peakCallingMACS2_annotationHomer.sh ${b} ${peaksdir} ${species} ${projName}; done;
cmd="parallel --tmpdir /media/rad/SSD1/atac_temp ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Run multiqc on mapped data and peaks
multiqc -o ${multiqcDir} -n ${projName}_stats ${projDir}


# 3) Christine Klement
# 3.1) For the first run for AGSaur_Snail001 sampels
jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="mm10"
user="christine"
projName="AGSaur_Snail001"
fastqdir="/media/rad/HDD1/atacseq/christine/AGSaur_Snail001/fastq"
outputDir="/media/rad/HDD1/atacseq/christine/AGSaur_Snail001"
# outputDir="/media/rad/SSD1/atac_temp/christine/AGSaur_Snail001"
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/${projName}"
multiqcDir="${outputDir}/qc/multiqc"; mkdir -p ${multiqcDir}
fastqcDir="${outputDir}/qc/fastqc"; mkdir -p ${fastqcDir}/00_original
logsDir="${outputDir}/logs"; mkdir -p ${logsDir}

bash scripts/01_map_singleendFastQ_bowtie2.sh ${fastqdir} ${outputDir} ${projName} ${species}
cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Run multiqc on original files
ls ${fastqdir}/*.gz | parallel --progress --eta -j 64 "fastqc -o ${fastqcDir}/00_original {}"
multiqc -o ${multiqcDir} -n 01_original_fastq_stats_${projName} ${fastqcDir}/00_original

# Run multiqc on trimmed files
multiqc -o ${multiqcDir} -n 02_trimmed_fastq_stats_${projName} ${fastqcDir}/01_trimmed ${logsDir}/trimmingLogs

# Get mapping stats and run multiqc on mapping logs
bamdir="${outputDir}/bams/trimmed"
ls ${bamdir}/*.bam | parallel --progress --eta -j 16 "samtools stats -@ 32 {} > {.}_mapping.logs"
mappingLogsDir="${outputDir}/logs/mappingLogs"; mkdir -p ${mappingLogsDir}; mv ${bamdir}/*.logs ${mappingLogsDir}
multiqc -o ${multiqcDir} -n 03_mapping_stats_${projName} ${mappingLogsDir}

# # Get peaks and homer annotation
# jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
# species="mm10"
# user="christine"
# projName="AGRad_ATACseq_MUC001"
# fastqdir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/demultiplexing/output/fastq"
# outputDir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001"
# bamdir="${outputDir}/bams/trimmed"
# peaksdir="${outputDir}/peaks"
# scriptsdir="${jobdir}/scripts/02_peakCallingMACs_annotationHomer/${projName}"
# mkdir -p ${scriptsdir}
# for b in ${bamdir}/*.bam; do bash scripts/02_peakCallingMACS2_annotationHomer.sh ${b} ${peaksdir} ${species} ${projName}; done;
# cmd="parallel --tmpdir /media/rad/SSD1/atac_temp ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}
