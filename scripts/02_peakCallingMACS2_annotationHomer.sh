#!/bin/bash

# USAGE: bash 02_peakCallingMACS2_annotationHomer.sh output/hgStomachF35/ENCFF300KHX_atacseq_hg_R1_PE100nt.bam output/hgStomachF35 hg19

# Get input arguments
jobdir="/home/rad/users/gaurav/projects/seqAnalysis/atacseq"
trtBam=${jobdir}/$1 # output/hgStomachF35/ENCFF300KHX_atacseq_hg_R1_PE100nt_subset1000.bam
outdir=${jobdir}/$2 # output/hgStomachF35
species=$3          # hg19 or mm10
projName=$4         # hgStomachF35

# Set user defined environment variables
bname=$(basename ${trtBam} .bam)
peaksoutdir=${outdir}/${bname}/macs2peaks
gooutdir=${outdir}/${bname}/GOterms
motifoutdir=${outdir}/${bname}/motifs
logdir=${outdir}/${bname}/logs
genomeFasta="${jobdir}/input/annotation/gtf/${species}_genes.gtf"
genomeFasta="${jobdir}/input/annotation/genomeData/${species}.fa"
scriptsdir=${jobdir}/scripts/02_peakCallingMACs_annotationHomer/${projName}
logFile="${logdir}/${bname}_peakCallingMacs2.log" 
homer="/home/rad/packages/homer/bin"
gSPC="human"
spc="hs"
if [ ${species} = "mm10" ]; then
 gSPC="mouse"
 spc="mm"
fi
if [ ${species} = "mm9" ]; then
 gSPC="mouse"
 spc="mm"
fi

# Create required dirs
mkdir -p ${scriptsdir} ${peaksoutdir} ${gooutdir} ${motifoutdir} ${logdir}

# Get the jobname to submit for each job

jobname="02_${bname}"

# Initialize a script file
scriptFile="${scriptsdir}/${bname}_peakCallingMacs2_homerAnnotation.sh"
touch ${scriptFile}
echo "#!/bin/bash" > "${scriptFile}"
echo "" >> "${scriptFile}"

# MACs Peak Calling
echo "# Perform MACs peak calling 2>&1 | tee ${logFile}" >> "${scriptFile}"
echo "macs2 callpeak --name ${bname} --treatment ${trtBam} --outdir ${peaksoutdir} --format BAM --pvalue 1e-5 --gsize ${spc} 2>&1 | tee -a ${logFile}" >> "${scriptFile}"
echo "" >> "${scriptFile}"

# Peak Annotation Using Homer
peaksFile=${peaksoutdir}/${bname}_peaks.narrowPeak
annPeaks=${peaksoutdir}/annotated_${bname}_peaks.txt
echo "# Perform peak annotation using Homer" >> "${scriptFile}"
echo "${homer}/annotatePeaks.pl ${peaksFile} ${species} > ${annPeaks} 2>&1 | tee -a ${logFile}" >> "${scriptFile}"
echo "" >> "${scriptFile}"

# GO Term Enrichment Using Homer
entrezidsFile=${peaksoutdir}/entrezIDs_${bname}.txt
echo "# Perform GO term enrichment using Homer 2>&1 | tee -a ${logFile}" >> "${scriptFile}"
echo "awk -F'\t' 'NR > 1 {print \$12}' ${annPeaks} > ${entrezidsFile} 2>&1 | tee -a ${logFile}" >> "${scriptFile}"
echo "${homer}/findGO.pl ${entrezidsFile} ${gSPC} ${gooutdir} 2>&1 | tee -a ${logFile}" >> "${scriptFile}"
echo "" >> "${scriptFile}"

# Motif Finding Using Homer
echo "# Perform motif finding using Homer" >> "${scriptFile}"
echo "${homer}/findMotifsGenome.pl ${annPeaks} ${species} ${motifoutdir} -len 6,10,15  2>&1 | tee -a ${logFile}" >> "${scriptFile}"
echo "" >> "${scriptFile}"

# queue="fat"
# chmod 775 "${scriptFile}"
# echo "Submitting job: ${jobname} ..." 
# if [[ $queue =~ "short" ]]
# then 
# qtime="01:59"
# else
# qtime="44:00"
# fi
# bsub -q ${queue} -n 6 -W ${qtime} -e "${errorsFile}" -o "${stdoutFile}" -J ${jobname} sh "${scriptFile}"

# Write the command in the script file and give it correct permission to run
chmod 775 "${scriptFile}"
echo "${scriptFile}"
echo " "



##################################### DOCS ########################################
#MACS Command Line
#An example of a macs command line with commonly used options is below
#macs2 callpeak  \
	#--name macs1
	#--treatment file1.bam file2.bam \
	#--control   ctrl1.bam ctrl2.bam  \
	#--outdir macs1 \
	#--format BAM \
	#--name macs1 \
	#--pvalue 1e-5 \
	#--gsize hs \
	#--tsize 75

#--treatment	The condition alignments in sam format.  Multiple files can be specified separated by spaces.
#--control	The control alignments in sam format
#--pvalue	Use 1e-5 at least
#--name	Experiment name, which will be used to generate output file names. DEFAULT: "NA"
#--outdir	The output directory where the results will go
#--gsize	Genome size.  Use hs for human or the actual base pair size
#--format	SAM (can also use BAMs and other formats)
#--tsize	Read length
#--bdg	Create a bed graph  (optional)
#--wig	Create a wig file  (optional)

# The only essential parameter is the --treatment one but I recommend you include the others (up to tsize) as a record of how you ran the program. More parameters and their descriptions are at https://github.com/taoliu/MACS/blob/macs_v1/README.rst

# Output files: The main output file is  NAME_peaks.narrowpeaks which is a  BED format file which contains the peak locations.   This will be the file you use as input for the motif analysis


# Peak Annotation Using Homer
#One of the many things that Homer can do is take a set of regions in a genome (our peaks for example) and annotate them with respect to where they lie compared to genes.  Each peak is assigned to a gene using the closest TSS (transcription start site).
#The command looks like
#annotatePeaks.pl \
#<peakfile> \
#<genome>   \
#> <outputfile>
#You can load the output file into excel if you want to look at your individual peaks later.
#Output columns .
#Description of Columns:

#1.  Peak ID
#2.  Chromosome
#3.  Peak start position
#4.  Peak end position
#5.  Strand
#6.  Peak Score
#7.  FDR/Peak Focus Ratio/Region Size
#8.  Annotation (i.e. Exon, Intron, ...)
#9.  Detailed Annotation (Exon, Intron etc. + CpG Islands, repeats, etc.)
#10.  Distance to nearest RefSeq TSS
#11.  Nearest TSS: Native ID of annotation file
#12.  Nearest TSS: Entrez Gene ID
#13.  Nearest TSS: Unigene ID
#14.  Nearest TSS: RefSeq ID
#15.  Nearest TSS: Ensembl ID
#16.  Nearest TSS: Gene Symbol
#17.  Nearest TSS: Gene Aliases
#18.  Nearest TSS: Gene description
#19.  Additional columns depend on options selected when running the program.

#Motif Finding using Homer
#http://homer.salk.edu/homer/introduction/basics.html Motif finding has many pitfalls.
#The motifs are small and common in the genome,
#The motifs are not always perfectly described in the database
#There is a lot of overlap beween motifs (degeneracy)
#Methods look for differential enrichment between a condition and control set of sequences.
#Controls are hard to get right and need thought.
#By default HOMER will use confident, non-regulated promoters as background when analyzing promoters, and sequences in the vicinity of genes for ChIP-Seq analysis (i.e. from –50kb to +50kb).  In each case sequences are matched for their GC content to avoid bias from CpG Islands.

#Homer motif finding usage http://homer.salk.edu/homer/ngs/peakMotifs.html

#findMotifsGenome.pl <mypeakfile>  \
#<mygeneome|mygenomefile>  \
#<outputdir> \
#-size <regionsize> \
#-p <threads> \
#-len <len1>,<len2>,<len3>

#- mypeakfile  	 the peak file from the macs output (.narrowpeak)
#- mygenome    	 hg19,mm9
#- mygenomefile	 hg19.fa mm9.fa (it will preprocess these)
#- outputdir   	 output directory
#- size        	 the size of the region (default is the whole peak)  (start with 50 for finding the main factor and work out)
#- p           	 number of threads
#- len         	 the lengths of the motifs to consider separated by commas
#- S           	 number of motifs to find.  The default is 25 which the docs say is pretty high already
#- dumpFasta   	 dumps the fasta of the regions

#########################################################
# 1 = PeakID (cmd=annotatePeaks.pl /usr/users/gjain/bin/projects/atacseq/output/hgStomachF35/peaks/ENCFF300KHX_atacseq_hg_R1_PE100nt_peaks.narrowPeak hg19)
# 2 = Chr
# 3 = Start
# 4 = End
# 5 = Strand
# 6 = Peak Score
# 7 = Focus Ratio/Region Size
# 8 = Annotation
# 9 = Detailed Annotation
#10 = Distance to TSS
#11 = Nearest PromoterID
#12 = Entrez ID
#13 = Nearest Unigene
#14 = Nearest Refseq
#15 = Nearest Ensembl
#16 = Gene Name
#17 = Gene Alias
#18 = Gene Description
#19 = Gene Type


