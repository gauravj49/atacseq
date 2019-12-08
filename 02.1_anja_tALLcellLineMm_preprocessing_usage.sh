# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Copy files to the HDD1
ls /media/nas/raw/PUB_CRCs/download_5/atacPE/*.gz | parallel --progress --eta -j 32 'rsync -arvP {} /media/rad/HDD1/atacseq/anja/tALLcellLineMm/fastq/pe'

# Rename the files with new naming convention complianced with illumina naming convention
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
