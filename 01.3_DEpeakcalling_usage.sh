# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Add Tobias Rausch's atacseq pipeline
cd scripts
git clone https://github.com/tobiasrausch/ATACseq.git 
cd ATACseq
sudo make all

outputDir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001"
# Differential peak calling
# Merge peaks across samples and create a raw count matrix.
DEDir="${outputDir}/de"; 
deInputDir=${DEDir}/input; mkdir -p ${deInputDir} 
deOutputDir=${DEDir}/output; mkdir -p ${deOutputDir} 
ls ${outputDir}/peaks/*/macs2peaks/*_summits.bed > ${deInputDir}/peaks_list.txt
ls ${outputDir}/bams/trimmed/*.bam > ${deInputDir}/bams_list.txt

# Remove samples
# 004_atac_* which are from a different run
# 005_atac_033	53646_LivMet-1

bash scripts/ATACseq/src/count.sh mm10 ${deInputDir}/peaks_list.txt ${deInputDir}/bams_list.txt AGRad_ATACseq_MUC001_peaks
mv *.gz ${deInputDir}



# To call differential peaks on a count matrix for TSS peaks, called counts.tss.gz, using DESeq2 we first need to create a file with sample level information (sample.info). For instance, if you have 2 replicates per condition:
# echo -e "name\tcondition" > sample.info
# zcat counts.tss.gz | head -n 1 | cut -f 5- | tr '\t' '\n' | sed 's/.final$//' | awk '{print $0"\t"int((NR-1)/2);}' >> sample.info
zcat ${deInputDir}/*tss.counts.gz | head -n 1 | cut -f 5- | tr '\t' '\n' | sed 's/.final$//' | awk '{print $0"\t"int((NR-1)/2);}' >> ${deInputDir}/sample.info

# Change all conditions to unique manually (no replicates)

Rscript scripts/ATACseq/R/dpeaks.R ${deInputDir}/*tss.counts.gz ${deInputDir}/sample.info


