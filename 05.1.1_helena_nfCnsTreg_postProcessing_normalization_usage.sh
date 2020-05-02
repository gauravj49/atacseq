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

# Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
consensusPeaksBed="${analysisDir}/${projName}_$(basename ${analysisDir})_consensus_peaks.bed"
rawCountsTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_rawCounts.txt"
peaksAnnTxtFile="${analysisDir}/$(basename ${consensusPeaksBed} .bed)_annotation.txt"

# 1) Parse concensus raw matrix and boolean matrix to get annoation files
echo "bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}"
bash scripts/parse_nfatac_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${outdir} ${origConsFile} ${origAnnFile} ${jobdir}

# Output files are:
echo "- Consensus bed file: ${consensusPeaksBed}"
echo "- Raw peaks count   : ${rawCountsTxtFile}"
echo "- Peaks annotation  : ${peaksAnnTxtFile}"

# 2) Normalize raw matrix with vst