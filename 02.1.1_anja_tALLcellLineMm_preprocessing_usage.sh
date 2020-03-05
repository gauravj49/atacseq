# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Copy files to the HDD1
ls /media/nas/raw/PUB_CRCs/download_5/atacPE/*.gz | parallel --progress --eta -j 32 'rsync -arvP {} /media/rad/HDD1/atacseq/anja/tALLcellLineMm/fastq/pe'

# Rename the files with new naming convention compliance with Illumina naming convention
# Naming convention: ExperimentInfo_Species_SeqProtocol_SeqType_Read
# 1) Experiment info: ExperimentName-CellLine-SpecialCondition-ReplicateNumber
# 2) Species: Human, mouse, etc.
# 3) SeqProtocol: Rnaseq, ChIPseq etc.
# 4) SeqType: paired end or single end
# 5) Read: R1 for side 1 and R2 for side 2
# Example: 
#   - TransPB-CD4-CD4PosTcells-Rep3_mm_atacseq_se_R1.fastq.gz
#   - TransPB-CD4-GSM2056292-Rep1_mm_atacseq_pe_R1.fastq.gz and TransPB-CD4-GSM2056292-Rep1_mm_atacseq_pe_R2.fastq.gz

# TransPB = Transposon PiggyBac
mv BMB1_GSM2056292_mm_atacseq_R1_1.fastq.gz TransPB-Bcell-BoneMarrow-Rep1_mm_atacseq_pe_R1.fastq.gz
mv BMB1_GSM2056292_mm_atacseq_R1_2.fastq.gz TransPB-Bcell-BoneMarrow-Rep1_mm_atacseq_pe_R2.fastq.gz
mv BMB2_GSM2056293_mm_atacseq_R1_1.fastq.gz TransPB-Bcell-BoneMarrow-Rep2_mm_atacseq_pe_R1.fastq.gz
mv BMB2_GSM2056293_mm_atacseq_R1_2.fastq.gz TransPB-Bcell-BoneMarrow-Rep2_mm_atacseq_pe_R2.fastq.gz
mv BMCLP1_GSM2056294_mm_atacseq_R1_1.fastq.gz TransPB-CLP-BoneMarrow-Rep1_mm_atacseq_pe_R1.fastq.gz
mv BMCLP1_GSM2056294_mm_atacseq_R1_2.fastq.gz TransPB-CLP-BoneMarrow-Rep1_mm_atacseq_pe_R2.fastq.gz
mv BMCLP2_GSM2056295_mm_atacseq_R1_1.fastq.gz TransPB-CLP-BoneMarrow-Rep2_mm_atacseq_pe_R1.fastq.gz
mv BMCLP2_GSM2056295_mm_atacseq_R1_2.fastq.gz TransPB-CLP-BoneMarrow-Rep2_mm_atacseq_pe_R2.fastq.gz
mv BMHSC1_GSM2056296_mm_atacseq_R1_1.fastq.gz TransPB-HSC-BoneMarrow-Rep1_mm_atacseq_pe_R1.fastq.gz
mv BMHSC1_GSM2056296_mm_atacseq_R1_2.fastq.gz TransPB-HSC-BoneMarrow-Rep1_mm_atacseq_pe_R2.fastq.gz
mv BMHSC2_GSM2056297_mm_atacseq_R1_1.fastq.gz TransPB-HSC-BoneMarrow-Rep2_mm_atacseq_pe_R1.fastq.gz
mv BMHSC2_GSM2056297_mm_atacseq_R1_2.fastq.gz TransPB-HSC-BoneMarrow-Rep2_mm_atacseq_pe_R2.fastq.gz
mv BMMPP1_GSM2056304_mm_atacseq_R1_1.fastq.gz TransPB-MPP-BoneMarrow-Rep1_mm_atacseq_pe_R1.fastq.gz
mv BMMPP1_GSM2056304_mm_atacseq_R1_2.fastq.gz TransPB-MPP-BoneMarrow-Rep1_mm_atacseq_pe_R2.fastq.gz
mv BMMPP2_GSM2056305_mm_atacseq_R1_1.fastq.gz TransPB-MPP-BoneMarrow-Rep2_mm_atacseq_pe_R1.fastq.gz
mv BMMPP2_GSM2056305_mm_atacseq_R1_2.fastq.gz TransPB-MPP-BoneMarrow-Rep2_mm_atacseq_pe_R2.fastq.gz
mv BMNK1_GSM2056308_mm_atacseq_R1_1.fastq.gz TransPB-NK-BoneMarrow-Rep1_mm_atacseq_pe_R1.fastq.gz
mv BMNK1_GSM2056308_mm_atacseq_R1_2.fastq.gz TransPB-NK-BoneMarrow-Rep1_mm_atacseq_pe_R2.fastq.gz
mv BMNK2_GSM2056309_mm_atacseq_R1_1.fastq.gz TransPB-NK-BoneMarrow-Rep2_mm_atacseq_pe_R1.fastq.gz
mv BMNK2_GSM2056309_mm_atacseq_R1_2.fastq.gz TransPB-NK-BoneMarrow-Rep2_mm_atacseq_pe_R2.fastq.gz
mv DN2aT1_GSM2692173_mm_atacseq_R1_1.fastq.gz TransPB-DN2a-Tcell-Rep1_mm_atacseq_pe_R1.fastq.gz
mv DN2aT1_GSM2692173_mm_atacseq_R1_2.fastq.gz TransPB-DN2a-Tcell-Rep1_mm_atacseq_pe_R2.fastq.gz
mv DN2aT2_GSM2692174_mm_atacseq_R1_1.fastq.gz TransPB-DN2a-Tcell-Rep2_mm_atacseq_pe_R1.fastq.gz
mv DN2aT2_GSM2692174_mm_atacseq_R1_2.fastq.gz TransPB-DN2a-Tcell-Rep2_mm_atacseq_pe_R2.fastq.gz
mv DN2bT1_GSM2692175_mm_atacseq_R1_1.fastq.gz TransPB-DN2b-Tcell-Rep1_mm_atacseq_pe_R1.fastq.gz
mv DN2bT1_GSM2692175_mm_atacseq_R1_2.fastq.gz TransPB-DN2b-Tcell-Rep1_mm_atacseq_pe_R2.fastq.gz
mv DN2bT2_GSM2692176_mm_atacseq_R1_1.fastq.gz TransPB-DN2b-Tcell-Rep2_mm_atacseq_pe_R1.fastq.gz
mv DN2bT2_GSM2692176_mm_atacseq_R1_2.fastq.gz TransPB-DN2b-Tcell-Rep2_mm_atacseq_pe_R2.fastq.gz
mv DN3T1_GSM2692177_mm_atacseq_R1_1.fastq.gz TransPB-DN3-Tcell-Rep1_mm_atacseq_pe_R1.fastq.gz
mv DN3T1_GSM2692177_mm_atacseq_R1_2.fastq.gz TransPB-DN3-Tcell-Rep1_mm_atacseq_pe_R2.fastq.gz
mv DN3T2_GSM2692333_mm_atacseq_R1_1.fastq.gz TransPB-DN3-Tcell-Rep2_mm_atacseq_pe_R1.fastq.gz
mv DN3T2_GSM2692333_mm_atacseq_R1_2.fastq.gz TransPB-DN3-Tcell-Rep2_mm_atacseq_pe_R2.fastq.gz
mv DN4T1_GSM2692178_mm_atacseq_R1_1.fastq.gz TransPB-DN4-Tcell-Rep1_mm_atacseq_pe_R1.fastq.gz
mv DN4T1_GSM2692178_mm_atacseq_R1_2.fastq.gz TransPB-DN4-Tcell-Rep1_mm_atacseq_pe_R2.fastq.gz
mv DN4T2_GSM2692179_mm_atacseq_R1_1.fastq.gz TransPB-DN4-Tcell-Rep2_mm_atacseq_pe_R1.fastq.gz
mv DN4T2_GSM2692179_mm_atacseq_R1_2.fastq.gz TransPB-DN4-Tcell-Rep2_mm_atacseq_pe_R2.fastq.gz
mv ETPT1_GSM2692171_mm_atacseq_R1_1.fastq.gz TransPB-ETP-Tcell-Rep1_mm_atacseq_pe_R1.fastq.gz
mv ETPT1_GSM2692171_mm_atacseq_R1_2.fastq.gz TransPB-ETP-Tcell-Rep1_mm_atacseq_pe_R2.fastq.gz
mv ETPT2_GSM2692172_mm_atacseq_R1_1.fastq.gz TransPB-ETP-Tcell-Rep2_mm_atacseq_pe_R1.fastq.gz
mv ETPT2_GSM2692172_mm_atacseq_R1_2.fastq.gz TransPB-ETP-Tcell-Rep2_mm_atacseq_pe_R2.fastq.gz
mv SPCD41_GSM2056332_mm_atacseq_R1_1.fastq.gz TransPB-CD4-SinglePositive-Rep1_mm_atacseq_pe_R1.fastq.gz
mv SPCD41_GSM2056332_mm_atacseq_R1_2.fastq.gz TransPB-CD4-SinglePositive-Rep1_mm_atacseq_pe_R2.fastq.gz
mv SPCD42_GSM2056333_mm_atacseq_R1_1.fastq.gz TransPB-CD4-SinglePositive-Rep2_mm_atacseq_pe_R1.fastq.gz
mv SPCD42_GSM2056333_mm_atacseq_R1_2.fastq.gz TransPB-CD4-SinglePositive-Rep2_mm_atacseq_pe_R2.fastq.gz
mv SPCD81_GSM2056306_mm_atacseq_R1_1.fastq.gz TransPB-CD8-BoneMarrow-Rep1_mm_atacseq_pe_R1.fastq.gz
mv SPCD81_GSM2056306_mm_atacseq_R1_2.fastq.gz TransPB-CD8-BoneMarrow-Rep1_mm_atacseq_pe_R2.fastq.gz
mv SPCD82_GSM2056307_mm_atacseq_R1_1.fastq.gz TransPB-CD8-SinglePositive-Rep2_mm_atacseq_pe_R1.fastq.gz
mv SPCD82_GSM2056307_mm_atacseq_R1_2.fastq.gz TransPB-CD8-SinglePositive-Rep2_mm_atacseq_pe_R2.fastq.gz

# Map the samples
jobdir='/home/rad/users/gaurav/projects/seqAnalysis/atacseq'
species="mm10"
user="anja"
projName="tALLcellLineMm"
projDir="/media/rad/HDD1/atacseq/anja/tALLcellLineMm"
fastqdir="${projDir}/fastq/pe"
bamdir="${projDir}/bams/trimmed"; mkdir -p ${bamdir}
bigwigdir="${bamdir}/bigwig"; mkdir -p ${bigwigdir}
scriptsdir="${jobdir}/scripts/01_map_fastQ_bowtie2/${projName}"
multiqcDir="${projDir}/qc/multiqc"; mkdir -p ${multiqcDir}
fastqcDir="${projDir}/qc/fastqc"; mkdir -p ${fastqcDir}
logsDir="${projDir}/logs"; mkdir -p ${logsDir}

bash scripts/01_map_pairedendFastQ_bowtie2.sh ${fastqdir} ${projDir} ${projName} ${species}
cmd="parallel ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Convert bam to bigwig using deeptools bamCoverage
# https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
effectiveGenomeSize=2652783500 # for mm10
for f in ${bamdir}/*_rmdup.bam; 
do
  echo ${f}; 
  bname=$(basename ${f} .bam); 
  # bamCoverage -b ${f} -o ${bigwigdir}/${bname}.bw -p 64; 
  bamCoverage -b ${f} -o ${bigwigdir}/${bname}_SeqDepthNorm.bw --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize ${effectiveGenomeSize} --ignoreForNormalization chrX --extendReads -p 64
  echo ""; 
done



# Get peaks and homer annotation
jobdir=" /home/rad/users/gaurav/projects/seqAnalysis/atacseq"
species="mm10"
user="anja"
projName="tALLcellLineMm"
projDir="/media/rad/HDD1/atacseq/anja/tALLcellLineMm"
fastqdir="${projDir}/fastq/pe"
mappingDir="${projDir}/mapping"
bamdir="${projDir}/bams/trimmed"
peaksdir="${projDir}/peaks"
scriptsdir="${jobdir}/scripts/02_peakCallingMACs_annotationHomer/${projName}"
mkdir -p ${scriptsdir} ${bam} ${peaksdir}
for b in ${bamdir}/*.bam; do bash scripts/02_peakCallingMACS2_annotationHomer.sh ${b} ${peaksdir} ${species} ${projName}; done;
cmd="parallel --tmpdir /media/rad/SSD1/atac_temp ::: "; for s in ${scriptsdir}/*.sh; do chmod 775 ${s}; cmd=$(echo "${cmd} ${s}"); done; eval ${cmd}

# Run multiqc on mapped data and peaks
multiqc -o ${multiqcDir} -n ${projName}_stats ${projDir}

# 4) Convert bam to bigwig using deeptools bamCoverage
bigwigdir="/media/rad/HDD1/temp_chip/snailchipP1/bigwig"; mkdir -p ${bigwigdir}
for f in /media/rad/HDD1/temp_chip/snailchipP1/mapping/*_rmdup.bam;
do
 echo ${f}
 bname=$(basename ${f} .bam)
 bamCoverage -b ${f} -o ${bigwigdir}/${bname}.bw -p 64
 echo ""
done



# ExperimentName=Meyloid
# CellLine=CMP
# SpecialCondition=GSM1441272
# ReplicateNumber=Rep1