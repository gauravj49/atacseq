# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Get merged peaks
# Get analysis dirs
mkdir -p /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis

# modify the peaks list and output file in the megrePeaks.sh file at the begining
# List of summit.bed is in array in the mergePeaks.sh script
bash scripts/mergePeaks.sh 

# Get the tab seaprated raw counts of the merged peaks for all the samples
# Using deeptools
multiBamSummary BED-file --BED /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks.bed --bamfiles /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/bams/trimmed/*.bam --smartLabels -out /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.npz --outRawCounts /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.tab -p 64

# Add peaknames to the file
ipython
pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

input_file  = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.tab"
output_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.matrix"
jmat = pd.read_csv(input_file, sep="\t")
# Fix column names: From #chr to chr and remove <'>
jmat.columns = jmat.columns.str.strip().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('#', '').str.replace("\'", '').str.replace("_r1_001_rmdup", '')
# rename the first column from #chr to chr
jmat.rename(columns={ jmat.columns[0]: "chr" }, inplace = True)
# Add peaks names to the dataframe
jmat.insert (3, "name", ["atacPeak_{0}".format(x) for x in jmat.index.tolist()])
jmat['peakID'] = jmat['chr'].str.cat(jmat['start'].apply(str), sep='_').str.cat(jmat['end'].apply(str), sep='_').str.cat(jmat['name'], sep='_')
# Get column names 
colNames = jmat.columns.tolist()
# Move peakID to the front
colNames.insert(0, colNames.pop(colNames.index('peakID')))
# Reorder columns using df.reindex() function
jmat = jmat.reindex(columns= colNames)
# Drop additional columns
jmat.drop(columns=['chr','start','end','name'], inplace=True)
jmat.to_csv(output_file, index=False, header=True, sep="\t", float_format='%.0f')

# Get the top 1% most varying genes

# Principal component analysis
# Source: http://quinlanlab.org/tutorials/bedtools/bedtools.html

# Get summits peaks
mkdir -p /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/summitBed
rsync -arvP  /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/*/macs2peaks/*_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/summitBed

# Get Jaccard statistic for all pairwise comparisons
cd /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/summitBed
parallel "bedtools jaccard -a {1} -b {2} \
         | awk 'NR>1' \
         | cut -f 3 \
         > {1}.{2}.jaccard" \
         ::: `ls *summits.bed` ::: `ls *summits.bed`

# Create a single file containing the pairwise Jaccard measurements from all tests
find . \
    | grep jaccard \
    | xargs grep "" \
    | sed -e s"/\.\///" \
    | perl -pi -e "s/.bed./.bed\t/" \
    | perl -pi -e "s/.jaccard:/\t/" \
    > /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt

# Get proper names
cat /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt \
| sed -e 's/_R1_001_rmdup_summits.bed//g' \
> tmp.txt && mv tmp.txt /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt

# Convert pairwise interaction to matrix
cat /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt | python /home/rad/users/gaurav/projects/seqAnalysis/atacseq/scripts/make_matrix.py > /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix

# Get labels for each dataset
cut -f 1 /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix > /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_labels.txt

# Start up R.
R

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(wesanderson) # devtools::install_github("karthik/wesanderson")
blues <- colorRampPalette(c('dark blue', 'light blue'))
greens <- colorRampPalette(c('dark green', 'light green'))
reds <- colorRampPalette(c('pink', 'dark red'))
 
x <- read.table('/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix')
labels <- read.table('/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_labels.txt')
ngroups <- length(unique(labels))
pca <- princomp(x)

nb.cols <- dim(labels)[1]
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
pdf('/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_pca.pdf')
qplot(pca$scores[,1], pca$scores[,2], color=factor(labels[,1]),     geom="point", size=1) + scale_color_manual(values = mycolors) 
# qplot(pca$scores[,1], pca$scores[,2], color=factor(labels[,1]),     geom="point", size=1) + scale_color_manual(values = c(blues(4), greens(5), reds(5))) 
# qplot(pca$scores[,1], pca$scores[,2], color=factor(labels[,1]),     geom="point", size=1) + scale_color_manual(values = c(blues(5), greens(4), reds(6), oranges(3))) 
# qplot(pca$scores[,1], pca$scores[,2], color=factor(labels[,1]),     geom="point", size=1) + scale_color_manual(values=wes_palette(n=18, name="GrandBudapest2"))
# qplot(pca$scores[,1], pca$scores[,2], color=factor(labels[,1]),     geom="point", size=1) + scale_color_viridis(discrete = TRUE, option = "D")+ scale_fill_viridis(discrete = TRUE)  
dev.off()

# Plot similarity heatmap
ipython
sns.set(font_scale=0.3)
jmat = pd.read_csv("/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix", sep="\t", index_col=0) 
g = sns.clustermap(jmat, vmax=0.02, cmap="ocean"); 
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 5)
plt.savefig("/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_heatmap.pdf")               

