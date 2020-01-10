# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Run CREAM on individual peaks
peaksDir='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/peaks'
outDir=$(dirname ${peaksDir})/cores/all; mkdir -p ${outDir}

for p in ${peaksDir}/*.narrowPeak;
do
  echo "Rscript scripts/R_identify_COREs_using_CREAM.R -if=${p} -of=${outDir}/$(basename ${p} .narrowPeak)_COREs.bed"
  Rscript scripts/R_identify_COREs_using_CREAM.R -if=${p} -of=${outDir}/$(basename ${p} .narrowPeak)_COREs.bed
  echo ""
done

# Get filtered cores with score > 10000
allcoresDir='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/all'
filcoresDir='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered'
for f in ${allcoresDir}/*.bed;
do
  bname=$(basename ${f} .bed); echo ${bname}
  oname=${filcoresDir}/${bname}_score_gt_10k.bed
  awk '{ if($5 >= 10000) { print }}' ${f} > ${oname}
done

# Sort the bed files
for f in *.bed; do echo ${f}; sort -k1,1 -k2n ${f} -o ${f}; done;

# Intersect the bed files for individual samples
filBedDir='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered'
cd ${filBedDir}
for m in 5320 53631 53646 6075;
do 
  multiIntersectBed -i ${m}*.bed -header | sort -nk4,4 -nk5,5 -nk1,1 | sed 's/_R1_001_rmdup_peaks_COREs_score_gt_10k.bed//g' >  common_cores_${m}.txt
done
cd -

# Get CRCs for the cores
cd /home/rad/users/gaurav/projects/ctrc/scripts/CLL_TFnetworks_2018
python2.7 CRC2.py -e /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered/unique_cores_5320_PPT-1_005.txt -b /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/bams/trimmed/5320_PPT-1_005_atac_028_S9_R1_001_rmdup.bam -g mm10 -o /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/crcs -n crcs_for_unique_cores_5320_PPT-1_005
cd -

# Get motifs using homer
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1 }' /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered/unique_cores_5320_PPT-1_005.txt > tmp && mv tmp /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered/unique_cores_5320_PPT-1_005.txt

awk '{ print $2 "\t" $3 "\t" $4 "\t" $1 }' /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered/all_common_cores_5320.txt > tmp && mv tmp /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered/all_common_cores_5320.txt

jobdir='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered'
species='mm10'

for n in  all_common_cores_5320 unique_cores_5320_PPT-1_005;
do 
  motifoutdir="/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/cores/filtered/motifs/${n}"
  f="${jobdir}/${n}.txt"
  findMotifsGenome.pl ${f} ${species} ${motifoutdir} -len 6,8,10  -p 8
done

