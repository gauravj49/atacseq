# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Parameters for the script
species="human"
user="helena"
projName="cnsTreg"
outdir="/media/rad/HDD1/atacseq"
origConsFile="/media/rad/HDD1/atacseq/helena/cnsTreg/results/bwa/mergedReplicate/macs/broadPeak/consensus/deseq2/consensus_peaks.mRp.clN.featureCounts.txt"
origAnnFile="/media/rad/HDD1/atacseq/helena/cnsTreg/results/bwa/mergedReplicate/macs/broadPeak/consensus/consensus_peaks.mRp.clN.boolean.txt"
jobdir="/home/rad/users/gaurav/projects/seqAnalysis/atacseq"

# Run the script
echo "bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}"
bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}
