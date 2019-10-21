#! /bin/bash
#
# STAR_HTSeq_arrays.batch
#
#SBATCH -J dfilterATAC
#SBATCH -p slim18
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -t 0-3:00
#SBATCH --mem 20000
#SBATCH -o findPeaks_%A.out
#SBATCH -e findPeaks_%A.err

#Finding peaks with Homer and merging the datasets

module load ngs/sratoolkit
module load ngs/samtools
module load ngs/Homer
module load ngs/bedtools2
module load R
module load ngs/UCSCutils

TDI=/work/project/becgsc_012/input_mice/Input.c_m1


TDlist=$(awk '{print $2}' < 01_samples_ATAC.txt | paste -s -d \ )

while read ID TD
do 

#Finding peaks
findPeaks $TD -style factor -o $ID\_ATAC.txt -i $TDI
sed 1,45d $ID\_ATAC.txt > $ID\_ATAC_2.txt
pos2bed.pl $ID\_ATAC_2.txt > $ID\_ATAC.bed

#filter out promoters
sortBed -i $ID\_ATAC.bed > $ID\_ATAC_sorted.bed

closestBed -a $ID\_ATAC_sorted.bed -b /work/project/becgsc_012/genome_mm10/Annotation/Genes/mm10_TSS_sorted.bed -d > $ID\_ATAC_mm10_TSS_distance.txt
Rscript 02_filter_promoters_Acinar.R $ID\_ATAC_mm10_TSS_distance.txt $ID\_ATAC_nonPromoter.bed

done < 01_samples_ATAC.txt

rm *_2.txt

#putting files together
cat *_nonPromoter.bed > multi_ATAC_nonPromoter.bed

#sorting
sortBed -i multi_ATAC_nonPromoter.bed > multi_ATAC_nonPromoter_sorted.bed

#merging files
bedtools merge -i multi_ATAC_nonPromoter_sorted.bed > multi_ATAC_nonPromoter_merged.bed

#annotatePeaks from Homer
annotatePeaks.pl multi_ATAC_nonPromoter_merged.bed mm10 -d $TDI $TDlist > 00_annotatePeaks_ATAC_nonPromoter.txt

#copy output files
mv findPeaks_$SLURM_JOB_ID.out 00_findPeaks_ATAC.out
mv findPeaks_$SLURM_JOB_ID.err 00_findPeaks_ATAC.err

