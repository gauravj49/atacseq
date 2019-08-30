# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

############### For Users ###################
# Get the wrappers and then run them in parallel

# 1) Christine Klement
# 1.1) For the first run for ATACseq_MUC001 sampels
species="mm10"
fastqdir="/home/rad/users/gaurav/projects/misc/output/AGRad_ATACseq_MUC001"
outputDir="output/AGRad_ATACseq_MUC001/mapping"
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/$(basename ${fastqdir})"

bash 01_map_singleendFastQ_bowtie2.sh ${fastqdir} ${outputDir} ${species}
cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Run multiqc
multiqcDir="${outputDir}/multiqc"; mkdir -p ${multiqcDir}
multiqc -o ${multiqcDir} -n AGRad_ATACseq_MUC001 ${outputDir}