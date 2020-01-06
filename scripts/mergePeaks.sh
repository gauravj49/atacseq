#!/bin/bash

# Taken from: https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/
# modified from https://bedops.readthedocs.io/en/latest/content/usage-examples/master-list.html

# summit_bed=(sample1_summits.bed sample2_summits.bed sample3_summits.bed)
# out=fAdrenal.master.merge.bed

# cd /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/summitBed
# ls $PWD/*.bed

summit_bed=(/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-Bcell-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-Bcell-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-CD4-SinglePositive-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-CD4-SinglePositive-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-CD8-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-CD8-SinglePositive-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-CLP-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-CLP-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN2a-Tcell-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN2a-Tcell-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN2b-Tcell-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN2b-Tcell-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN3-Tcell-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN3-Tcell-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN4-Tcell-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DN4-Tcell-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DP-Tcell-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-DP-Tcell-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-ETP-Tcell-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-ETP-Tcell-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-HSC-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-HSC-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-MPP-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-MPP-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-NK-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/broadPeaks/TransPB-NK-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup_peaks.broadPeak)

out=/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks.bed

tmpd=/tmp/tmp$$
mkdir -p $tmpd

## First, union all the peaks together into a single file.
bedlist=""
for bed in ${summit_bed[*]}
do
    bedlist="${bedlist} ${bed}"
done

# extended summits to 500-bp windows (±250 bp)
bedops --range 250 -u $bedlist > $tmpd/tmp.bed

## The master list is constructed iteratively.  For each pass through
## the loop, elements not yet in the master list are merged into
## non-overlapping intervals that span the union (this is just bedops
## -m).  Then for each merged interval, an original element of highest
## score within the interval is selected to go in the master list.
## Anything that overlaps the selected element is thrown out, and the
## process then repeats.
iters=1
solns=""
stop=0
while [ $stop == 0 ]
do
    echo "merge steps..."

    ## Condense the union into merged intervals. This klugey bit
    ## before and after the merging is because we don't want to merge
    ## regions that are simply adjacent but not overlapping
    bedops -m --range 0:-1 $tmpd/tmp.bed \
        | bedops -u --range 0:1 - \
        > $tmpd/tmpm.bed

    ## Grab the element with the highest score among all elements forming each interval.
    ## If multiple elements tie for the highest score, just grab one of them.
    ## Result is the current master list.  Probably don't need to sort, but do it anyway
    ## to be safe since we're not using --echo with bedmap call.
    bedmap --max-element $tmpd/tmpm.bed $tmpd/tmp.bed \
        | sort-bed - \
        > $tmpd/$iters.bed
    solns="$solns $tmpd/$iters.bed"
    echo "Adding `awk 'END { print NR }' $tmpd/$iters.bed` elements"

    ## Are there any elements that don't overlap the current master
    ## list?  If so, add those in, and repeat.  If not, we're done.
    bedops -n 1 $tmpd/tmp.bed $tmpd/$iters.bed \
       > $tmpd/tmp2.bed

    mv $tmpd/tmp2.bed $tmpd/tmp.bed

    if [ ! -s $tmpd/tmp.bed ]
    then
        stop=1
    fi

    ((iters++))
done

## final solution
bedops -u $solns \
   > $out

## Clean up
rm -r $tmpd

exit 0

















# summit_bed=(/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/*/macs2peaks/*_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup/macs2peaks/5320_LivMet-1_005_atac_030_S11_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_LivMet-3_005_atac_031_S12_R1_001_rmdup/macs2peaks/5320_LivMet-3_005_atac_031_S12_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_LungMet-1_005_atac_029_S10_R1_001_rmdup/macs2peaks/5320_LungMet-1_005_atac_029_S10_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_PPT-1_004_atac_023_20k_S18_R1_001_rmdup/macs2peaks/5320_PPT-1_004_atac_023_20k_S18_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/5320_PPT-1_005_atac_028_S9_R1_001_rmdup/macs2peaks/5320_PPT-1_005_atac_028_S9_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_LivMet-1_005_atac_025_S6_R1_001_rmdup/macs2peaks/53631_LivMet-1_005_atac_025_S6_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_LungMet-2_005_atac_026_S7_R1_001_rmdup/macs2peaks/53631_LungMet-2_005_atac_026_S7_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_LungMet-3_005_atac_027_S8_R1_001_rmdup/macs2peaks/53631_LungMet-3_005_atac_027_S8_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53631_PPT-1_005_atac_024_S5_R1_001_rmdup/macs2peaks/53631_PPT-1_005_atac_024_S5_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_LivMet-1_005_atac_033_S2_R1_001_rmdup/macs2peaks/53646_LivMet-1_005_atac_033_S2_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_LivMet-2_005_atac_034_S3_R1_001_rmdup/macs2peaks/53646_LivMet-2_005_atac_034_S3_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_LivMet-3_005_atac_035_S4_R1_001_rmdup/macs2peaks/53646_LivMet-3_005_atac_035_S4_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_PPT-1_004_atac_022_S17_R1_001_rmdup/macs2peaks/53646_PPT-1_004_atac_022_S17_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_PPT-1_005_atac_032_S1_R1_001_rmdup/macs2peaks/53646_PPT-1_005_atac_032_S1_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/53646_PPT-1_frozen_004_atac_022_fr_S16_R1_001_rmdup/macs2peaks/53646_PPT-1_frozen_004_atac_022_fr_S16_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/6075_blood_005_atac_038_S15_R1_001_rmdup/macs2peaks/6075_blood_005_atac_038_S15_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/6075_LungMet-1_005_atac_037_S14_R1_001_rmdup/macs2peaks/6075_LungMet-1_005_atac_037_S14_R1_001_rmdup_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/6075_PPT-1_005_atac_036_S13_R1_001_rmdup/macs2peaks/6075_PPT-1_005_atac_036_S13_R1_001_rmdup_summits.bed)

# summit_bed=(/media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS244_12h_NKT_TAAGGC_atacseq_mm_se_rmdup/macs2peaks/GS244_12h_NKT_TAAGGC_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS245_12h_NKT_CGTACT_atacseq_mm_se_rmdup/macs2peaks/GS245_12h_NKT_CGTACT_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS246_5d_NKT_AGGCAG_atacseq_mm_se_rmdup/macs2peaks/GS246_5d_NKT_AGGCAG_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS247_5d_NKT_TCCTGA_atacseq_mm_se_rmdup/macs2peaks/GS247_5d_NKT_TCCTGA_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS248_NKT1_GGACTC_atacseq_mm_se_rmdup/macs2peaks/GS248_NKT1_GGACTC_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS249_NKT1_TAGGCA_atacseq_mm_se_rmdup/macs2peaks/GS249_NKT1_TAGGCA_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS250_DP_CTCTCT_atacseq_mm_se_rmdup/macs2peaks/GS250_DP_CTCTCT_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS251_DP_CAGAGA_atacseq_mm_se_rmdup/macs2peaks/GS251_DP_CAGAGA_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS394_24h_NKT_TAAGGCGA_atacseq_mm_se_rmdup/macs2peaks/GS394_24h_NKT_TAAGGCGA_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS395_24h_NKT_TAAGGCGA_atacseq_mm_se_rmdup/macs2peaks/GS395_24h_NKT_TAAGGCGA_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS396_36h_NKT_CGTACTAG_atacseq_mm_se_rmdup/macs2peaks/GS396_36h_NKT_CGTACTAG_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS397_36h_NKT_CGTACTAG_atacseq_mm_se_rmdup/macs2peaks/GS397_36h_NKT_CGTACTAG_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS398_48h_NKT_AGGCAGAA_atacseq_mm_se_rmdup/macs2peaks/GS398_48h_NKT_AGGCAGAA_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS399_60h_NKT_TCCTGAGC_atacseq_mm_se_rmdup/macs2peaks/GS399_60h_NKT_TCCTGAGC_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS400_71h_NKT_GGACTCCT_atacseq_mm_se_rmdup/macs2peaks/GS400_71h_NKT_GGACTCCT_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS401_84h_NKT_TCCTGAGC_atacseq_mm_se_rmdup/macs2peaks/GS401_84h_NKT_TCCTGAGC_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS402_96h_NKT_TAGGCATG_atacseq_mm_se_rmdup/macs2peaks/GS402_96h_NKT_TAGGCATG_atacseq_mm_se_rmdup_peaks.broadPeak /media/rad/HDD1/atacseq/sabrina/nkTimecourse/peaks/GS403_NKT1_CTCTCTAC_atacseq_mm_se_rmdup/macs2peaks/GS403_NKT1_CTCTCTAC_atacseq_mm_se_rmdup_peaks.broadPeak)