# Principal component analysis
# Source: http://quinlanlab.org/tutorials/bedtools/bedtools.html

# Get Jaccard statistic for all pairwise comparisons
cd output/christine/AGRad_ATACseq_MUC001/analysis/summitBed
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
    > /home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt

# Get proper names
cat /home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt \
| sed -e 's/_R1_001_summits.bed//g' \
> tmp.txt && mv tmp.txt /home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt

# Convert pairwise interaction to matrix
cat /home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons.txt | python scripts/make_matrix.py > /home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix

# Get labels for each dataset
cut -f 1 /home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix > /home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_labels.txt

# Start up R.
R

library(ggplot2)
library(RColorBrewer)
blues <- colorRampPalette(c('dark blue', 'light blue'))
greens <- colorRampPalette(c('dark green', 'light green'))
reds <- colorRampPalette(c('pink', 'dark red'))
 
x <- read.table('/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix')
labels <- read.table('/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_labels.txt')
ngroups <- length(unique(labels))
pca <- princomp(x)

pdf('/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_pca.pdf')
# qplot(pca$scores[,1], pca$scores[,2], color=factor(labels[,1]),     geom="point", size=1) + scale_color_manual(values = c(blues(4), greens(5), reds(5))) 
qplot(pca$scores[,1], pca$scores[,2], color=factor(labels[,1]),     geom="point", size=1) + scale_color_manual(values = c(blues(5), greens(4), reds(6), blacks(3))) 
dev.off()



library(gplots)
library(RColorBrewer)
library(pheatmap)
jaccard_table <- x[, -1]
jaccard_matrix <- as.matrix(jaccard_table)

# pdf('/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_heatmap.pdf')
# heatmap.2(jaccard_matrix, col = brewer.pal(9,"Blues"), margins = c(14, 14), density.info = "none", lhei=c(2, 8), trace= "none")
corrHeatmapFile = '/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_heatmap.pdf'
pheatmap(jaccard_matrix, filename=corrHeatmapFile)
dev.off()

ipython
jmat = pd.read_csv("/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix", sep="\t", index_col=0) 
sns.clustermap(jmat, vmax=0.02, cmap="ocean");                                                                                                                                 plt.savefig("/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_heatmap.pdf")               


#  https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html

R
#

sortedBAM <- "/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/mapping/5320_LivMet-1_S11_R1_001.bam"


library(Rsubread)
# Takes ~ 15 minutes on 3.1 GHz Intel Core i7 Mac pro
pmapped <- propmapped(sortedBAM)
pmapped

library(Rsamtools)
library(ggplot2)
library(magrittr)
sortedBAMpng <- "/home/rad/users/gaurav/projects/seqAnalysis/atacseq/output/christine/AGRad_ATACseq_MUC001/analysis/5320_LivMet-1_S11_R1_001_mappedReadsDistribution.pdf"
pdf(sortedBAMpng);
idxstatsBam(sortedBAM) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + geom_bar(stat = "identity") + coord_flip()
dev.off()


library(ChIPseeker)

MacsCalls_chr20_filteredAnno <- annotatePeak(MacsCalls_chr20_filtered, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)