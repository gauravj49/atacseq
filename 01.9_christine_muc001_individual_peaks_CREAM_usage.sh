# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Run CREAM on individual peaks
peaksDir='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/individual/peaks'
outDir=$(dirname ${peaksDir})/cores; mkdir -p ${outDir}

for p in ${peaksDir}/*.narrowPeak;
do
  echo "Rscript scripts/R_identify_COREs_using_CREAM.R -if=${p} -of=${outDir}/$(basename ${p} .narrowPeak)_COREs.bed"
  Rscript scripts/R_identify_COREs_using_CREAM.R -if=${p} -of=${outDir}/$(basename ${p} .narrowPeak)_COREs.bed
  echo ""
done

