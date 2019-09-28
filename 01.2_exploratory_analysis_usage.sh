# Principal component analysis
# Source: http://quinlanlab.org/tutorials/bedtools/bedtools.html

# Get merged peaks
# modify the peaks list and output file in the megrePeaks.sh file at the begining
bash scripts/mergePeaks.sh 

# Get analysis dirs
mkdir -p /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis

# Get summits peaks
mkdir -p /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/summitBed
rsync -avrP rsync -arvP  /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/peaks/*/macs2peaks/*_summits.bed /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/summitBed

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
jmat = pd.read_csv("/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix", sep="\t", index_col=0) 
sns.clustermap(jmat, vmax=0.02, cmap="ocean");                                                                                                                                 plt.savefig("/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_heatmap.pdf")               

# Get the top 1% most varying genes


















#  https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html

R
#

sortedBAM <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/mapping/5320_LivMet-1_S11_R1_001.bam"


library(Rsubread)
# Takes ~ 15 minutes on 3.1 GHz Intel Core i7 Mac pro
pmapped <- propmapped(sortedBAM)
pmapped

library(Rsamtools)
library(ggplot2)
library(magrittr)
sortedBAMpng <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_LivMet-1_S11_R1_001_mappedReadsDistribution.pdf"
pdf(sortedBAMpng);
idxstatsBam(sortedBAM) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") + coord_flip()
dev.off()


library(ChIPseeker)

MacsCalls_chr20_filteredAnno <- annotatePeak(MacsCalls_chr20_filtered, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)