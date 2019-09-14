# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Create peaks dir
peaksdir="output/christine/AGRad_ATACseq_MUC001/peaks"
mkdir -p ${peaksdir}

# Make genome information (chromosome sizes) from the BAM file to get a genome.info file:
samtools view -H output/christine/AGRad_ATACseq_MUC001/mapping/5320_LivMet-1_S11_R1_001.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > ${peaksdir}/genome.info

# Run HMMRATAC on the sorted BAM ATACseq.sorted.bam, the BAM index file ATACseq.sorted.bam.bai, and the genome information file genome.info:
java -jar /home/rad/packages/hmmratac/HMMRATAC_V1.2.7_exe.jar -b output/christine/AGRad_ATACseq_MUC001/mapping/5320_LivMet-1_S11_R1_001.bam -i output/christine/AGRad_ATACseq_MUC001/mapping/5320_LivMet-1_S11_R1_001.bam.bai -g ${peaksdir}/genome.info



