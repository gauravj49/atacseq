# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Fastq Location: /media/rad/HDD1/PUB_CRCs/atac/atac_tALLcellLine_hg
# /home/rad/packages/data/fasta/human/hg19/hg19.fa

# Run analysis in parallel
cd output/anja/tALLcellLine_hg
for n in CCRFCEM_GSM3693106 MOLT4_GSM3693104 Jurkat_GSM3693103 RPMI8402_GSM3693105; do echo "/media/rad/HDD1/PUB_CRCs/atac/atac_tALLcellLine_hg/${n}_hs_atacseq_R1"; done | parallel --progress --eta -j 16 ' bash /home/rad/users/gaurav/projects/seqAnalysis/atacseq/scripts/ATACseq/src/atac.sh hg19 {}_1.fastq.gz {}_2.fastq.gz /home/rad/packages/data/fasta/human/hg19/hg19.fa {}'


bash /home/rad/users/gaurav/projects/seqAnalysis/atacseq/scripts/ATACseq/src/atac.sh hg19 /media/rad/HDD1/PUB_CRCs/atac/atac_tALLcellLine_hg/Jurkat_GSM3693103_hs_atacseq_R1_1.fastq.gz /media/rad/HDD1/PUB_CRCs/atac/atac_tALLcellLine_hg/Jurkat_GSM3693103_hs_atacseq_R1_2.fastq.gz /home/rad/packages/data/fasta/human/hg19/hg19.fa Jurkat_GSM3693103
