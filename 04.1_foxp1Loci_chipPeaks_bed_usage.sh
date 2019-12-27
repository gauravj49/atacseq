# pwd
cd /home/rad/users/gaurav/users/thorsten/foxp1Loci

# Get input foxp1 chipseq peaks files
# publication: Konopacki, C., Pritykin, Y., Rubtsov, Y., Leslie, C. S., & Rudensky, A. Y. (2019). Transcription factor Foxp1 regulates Foxp3 chromatin binding and coordinates regulatory T cell function. Nature immunology, 20(2), 232.
# Source: https://www.nature.com/articles/s41590-018-0291-z?platform=oscar&draft=collection#Sec30
# File: 41590_2018_291_MOESM3_ESM-2.xlsx
# Converted the treg excel to bed
# head /home/rad/users/gaurav/users/thorsten/foxp1Loci/input/foxp1_treg_peaks.bed 
# ┌───────┬──────────┬──────────┬──────────────────────────────────┐
# │ chrom │  start   │   end    │               name               │
# ├───────┼──────────┼──────────┼──────────────────────────────────┤
# │ chr17 │ 33716102 │ 33716770 │ foxp1_peak_24454_promoter_March2 │
# │ chr17 │ 33718550 │ 33719080 │ foxp1_peak_24456_promoter_March2 │
# │ chr15 │ 31530767 │ 31531403 │ foxp1_peak_19171_promoter_March6 │
# │ chr2  │ 60210280 │ 60210878 │ foxp1_peak_30703_promoter_March7 │
# └───────┴──────────┴──────────┴──────────────────────────────────┘

# 1) ################################################################################################################
ipython

# Get the peaks tab file for atacseq data
# Input and output files
input_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/5320_53631_53646_6075_merge_master_peaks_vstCounts.matrix"
peaksDF    = pd.read_csv(input_file, sep="\t", index_col=0)
# Save the topPeaks files in tab separated files
tabPeaksDF = peaksDF.copy()

# Get the chr, start and end columns
tabPeaksDF['tabDFSlice'] = tabPeaksDF.index.values.tolist()
tabDFSlice               = tabPeaksDF['tabDFSlice'].str.split("_", n = 3, expand = True)
tabPeaksDF['chr']        = tabDFSlice[0]
tabPeaksDF['start']      = tabDFSlice[1]
tabPeaksDF['end']        = tabDFSlice[2]

# Rearrange columns to new order
new_columns_order = ['chr', 'start', 'end', '5320_livmet-1_005_atac_030_s11', '5320_livmet-3_005_atac_031_s12', '5320_lungmet-1_005_atac_029_s10', '5320_ppt-1_004_atac_023_20k_s18', '5320_ppt-1_005_atac_028_s9', '53631_livmet-1_005_atac_025_s6', '53631_lungmet-2_005_atac_026_s7', '53631_lungmet-3_005_atac_027_s8', '53631_ppt-1_005_atac_024_s5', '53646_livmet-1_005_atac_033_s2', '53646_livmet-2_005_atac_034_s3', '53646_livmet-3_005_atac_035_s4', '53646_ppt-1_004_atac_022_s17', '53646_ppt-1_005_atac_032_s1', '53646_ppt-1_frozen_004_atac_022_fr_s16', '6075_blood_005_atac_038_s15', '6075_lungmet-1_005_atac_037_s14', '6075_ppt-1_005_atac_036_s13']
tabPeaksDF        = tabPeaksDF[new_columns_order]

# Save the DF as tab separated file
taboutput_file = "/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/5320_53631_53646_6075_merge_master_peaks_vstCounts.tab"
tabPeaksDF.to_csv(taboutput_file, sep='\t', index = False, float_format='%.2g')

############ USER DEFINIED CODE ############
def read_fasta(fp):
  ''' 
    - Parse a fasta file using a generator 
    - Source: https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
  '''
  name, seq = None, []
  for line in fp:
      line = line.rstrip()
      if line.startswith(">"):
          if name: yield (name, ''.join(seq))
          name, seq = line, []
      else:
          seq.append(line)
  if name: yield (name, ''.join(seq))

# 2) ################################################################################################################
# Intersect with atacseq peaks

# 1) Get the peaks in the foxp1 region
peakMatrix="/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/5320_53631_53646_6075_merge_master_peaks_vstCounts.tab"

# 2.2) Get bed files for each sample
od="$(dirname ${peakMatrix})/samplesBed"
mkdir -p ${od}
header=$(head -1 ${peakMatrix})
for i in {4..21};
do
 bname=$(cut -f${i} <<< ${header})
 echo ${bname}
 cut -f1-3,${i} ${peakMatrix} > ${od}/${bname}_peaks.bed
 # Sort the bed file 
 sort -k1,1 -k2n ${od}/${bname}_peaks.bed -o ${od}/${bname}_peaks.bed
done

# Remove random chromosome GL
for f in ${od}/*.bed; do sed -i '/GL/d' ${f}; done;
for f in ${od}/*.bed; do sed -i '/JH/d' ${f}; done;

# 2.3) Get the mean score of atacseq peaks for the foxp1 region
pd="$(dirname ${peakMatrix})/foxp1Peaks"; mkdir -p ${pd}
ib="/home/rad/users/gaurav/users/thorsten/foxp1Loci/input/foxp1_treg_peaks.bed"
sort -k1,1 -k2n ${ib} -o ${ib}
for bbed in ${od}/*.bed;
do
 echo ${bbed}
 bedtools map -a ${ib} -b ${bbed} -c 4 -o mean > ${pd}/$(basename ${bbed} .bed)_foxp1_peaks_mean.bed
 # Add header to the file
 header="chr\tstart\tend\tname\t$(basename ${bbed} .bed)"
 sed -i "1i${header}" ${pd}/$(basename ${bbed} .bed)_foxp1_peaks_mean.bed
done

# Merge the individual mean peaks into a dataframe
ipython
from scipy.stats import zscore

# Read and merge into df
pkDF = pd.concat([pd.read_csv(f, sep='\t').set_index(['chr', 'start', 'end', 'name']) for f in glob.glob('/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/foxp1Peaks/*.bed')],axis=1).reset_index()
pkDF.replace('.',0, inplace=True)

# Plot the heatmap
peaksDF = pkDF.copy()
peaksDF.drop(columns = ['chr','start','end'], inplace = True)
peaksDF.set_index('name',inplace=True)
peaksDF = peaksDF.loc[~(peaksDF==0).all(axis=1)] # Remove rows with all 0s
peaksDF = peaksDF.astype('float') # Conver all columns to float
# Remove _peaks.bed_mean from column names
peaksDF = peaksDF.rename(columns={col: col.split('_peaks')[0] for col in peaksDF.columns})

# Drop the 004 samples
peaksDF.drop(list(peaksDF.filter(regex = '004')), axis = 1, inplace = True)

# Save to
output_file = '/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/foxp1loci_merge_peaks_averageReads.txt'
peaksDF.to_csv(output_file, sep='\t', index = True, float_format='%.2g')

# Calculate the zscore of the dataframe
# The peaksDF.apply(zscore, axis=1) returns an array 
# https://stackoverflow.com/questions/35491274/pandas-split-column-of-lists-into-multiple-columns
from scipy.stats import zscore
peaksDF = peaksDF.astype('float')
zscorePeaksDF     = pd.DataFrame(pd.DataFrame(peaksDF.apply(zscore, axis=1))[0].values.tolist(), columns=peaksDF.columns.tolist(), index=peaksDF.index)

# Plot the heatmap for all replicates separately for all time points
sns.set(font_scale=1.25)
g = sns.clustermap(zscorePeaksDF, cmap='RdBu_r', figsize=(10,25), yticklabels=False); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0); cax = plt.gcf().axes[-1]; cax.tick_params(labelsize=5);

heatmapPlotPdf = "{0}_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')

#########################################
# Save session
import dill
filename = "{0}_session.pkl".format(get_file_info(output_file)[3])
dill.dump_session(filename)

# and to load the session again:
import dill
filename = "{0}_session.pkl".format(get_file_info(output_file)[3])
dill.load_session(filename)
#########################################

# Cluster information
# https://support.bioconductor.org/p/93424/
# https://stackoverflow.com/questions/27820158/pheatmap-in-r-how-to-get-clusters

# Sort input file and remove duplicate lines
f='/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/foxp1loci_merge_peaks_averageReads.txt'
(head -n 1 ${f} && tail -n +2 ${f} | sort -u) > tmp && mv tmp ${f}

# Draw heatmap with clusterIDs and save the clusterIDs in a separate file
R
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(data.table))

## Get the input data
cat("- Reading input file ...\n")
inputfile              <- '/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/foxp1loci_merge_peaks_averageReads.txt'
allDataOrig            <- data.frame(fread(inputfile, header=TRUE, sep="\t"), check.names=F)
row.names(allDataOrig) <- allDataOrig$name
allDataOrig[1]         <- NULL

# Get the number of rows and columns
ncols  <- length(names(allDataOrig))
nrows  <- length(rownames(allDataOrig))
cat("- Total samples : ", ncols,"\n- Total features: ", nrows, "\n")

# Create the output directories and files
cat("\n2) Create the output directories and files ...\n")
basefilename    <- tools::file_path_sans_ext(basename(inputfile))
extfilename     <- tools::file_ext(inputfile)
outfilebasename <- paste0(dirname(inputfile),"/",basefilename)

# Plot the heatmap with the new cluster and order and save it to the output file
pdffile         <- paste0(outfilebasename,"_clustered_peaks.pdf");
clusterfile     <- paste0(outfilebasename,"_mit_cluster_labels.txt");

# Get the clusters and save it to the file
numClusters          <- 4
# scaled_df            <- as.data.frame(matrix(scale(c(as.matrix(allDataOrig))), nrow=nrows, byrow=TRUE))
scaled_df            <- (allDataOrig - rowMeans(allDataOrig))/(rowSds(as.matrix(allDataOrig)))[row(allDataOrig)]
colnames(scaled_df)  <- colnames(allDataOrig)
row.names(scaled_df) <- row.names(allDataOrig)
resHeatmap           <- pheatmap(scaled_df, cutree_row=numClusters, show_rownames=F, silent=T)
allDataClust         <- cbind(scaled_df, cluster = cutree(resHeatmap$tree_row, k= numClusters))

# Save the clusterIDs to output file
fwrite(allDataClust, file=clusterfile, quote=F, row.names=T, sep='\t',  nThread=48)

# Get clusters for annoataion
annClust        <- allDataClust['cluster']

# Treat them as factors so that they are not plotted as continuous data
annClust$cluster <- as.factor(annClust$cluster)

# Draw the final cluster
finHeatmap      <- pheatmap(scaled_df, filename=pdffile, cutree_row=numClusters, show_rownames=F, annotation_row=annClust)
