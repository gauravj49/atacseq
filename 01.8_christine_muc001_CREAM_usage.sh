# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# Run CREAM
suppressPackageStartupMessages(library(CREAM))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(data.table))

IdentifiedCOREs <- CREAM(in_path = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks.bed", MinLength = 1000, peakNumMin = 2, WScutoff = 0.25 )

# > str(IdentifiedCOREs)
#  chr [1:1505, 1:6] "chr1" "chr1" "chr1" "chr1" "chr1" "chr1" "chr1" "chr1" ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : NULL
#   ..$ : NULL
# > 

# Save the data frame into a output file
outputFile       <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs.bed"
coreDF           <- data.frame(IdentifiedCOREs[,c(1:3,6)])
colnames(coreDF) <- c('chr','start','end','COREs_score')
fwrite(coreDF, file=outputFile, quote=F, row.names=F, col.names=F, sep='\t',  nThread=48)

# Sort file
system(paste0("sort -k1,1 -k2n ", outputFile, " -o ", outputFile))
# Cltr+D+D
#****************************************************************************************************

# Get the tab seaprated raw counts of the merged peaks for all the samples
# Using deeptools
multiBamSummary BED-file --BED "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs.bed" --bamfiles /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/bams/trimmed/*.bam --smartLabels -out "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_rawCounts.npz" --outRawCounts "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_rawCounts.tab" -p 64

#****************************************************************************************************
# Add peaknames to the file
ipython

pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

input_file  = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_rawCounts.tab"
output_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_rawCounts.matrix"
peaksDF = pd.read_csv(input_file, sep="\t")
# Fix column names: From #chr to chr and remove <'>
peaksDF.columns = peaksDF.columns.str.strip().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('#', '').str.replace("\'", '').str.replace("_r1_001_rmdup", '')
# rename the first column from #chr to chr
peaksDF.rename(columns={ peaksDF.columns[0]: "chr" }, inplace = True)
# Add peaks names to the dataframe
peaksDF.insert (3, "name", ["atacPeak_{0}".format(x) for x in peaksDF.index.tolist()])
peaksDF['peakID'] = peaksDF['chr'].str.cat(peaksDF['start'].apply(str), sep='_').str.cat(peaksDF['end'].apply(str), sep='_').str.cat(peaksDF['name'], sep='_')
# Get column names 
colNames = peaksDF.columns.tolist()
# Move peakID to the front
colNames.insert(0, colNames.pop(colNames.index('peakID')))
# Reorder columns using df.reindex() function
peaksDF = peaksDF.reindex(columns= colNames)
# Drop additional columns
peaksDF.drop(columns=['chr','start','end','name'], inplace=True)
peaksDF.to_csv(output_file, index=False, header=True, sep="\t", float_format='%.0f')
# Cltr+D+D
#****************************************************************************************************

# Normalize raw matrix with vst
R
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
inputFile  <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_rawCounts.matrix"
outputFile <- "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_vstCounts.matrix"
# Import data from featureCounts
countdata <- read.table(inputFile, header=TRUE, row.names=1, check.names = FALSE)
coldata   <- data.frame(row.names=colnames(countdata), samples=colnames(countdata))
dds       <- DESeqDataSetFromMatrix(countData=countdata, design=~1, colData=coldata)
vst         <- varianceStabilizingTransformation(dds, blind=TRUE)  
vsd         <- assay(vst)
# Convert rownames as feature column to save it in the csv file
vsd           <- cbind(feature = rownames(vsd), vsd)
rownames(vsd) <- 1:nrow(vsd)
# Save the counts file to output csv file
write.table(vsd   , file = outputFile , row.names = F, sep = '\t', quote = F)
# Cltr+D+D
#****************************************************************************************************

# Plot the pca for COREs identified by CREAM algorithm
ipython

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

# Input and output files
input_file      = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_vstCounts.matrix"
rankoutput_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_vstCounts.ranks"

# Import data into the dataframe
originalPeaksDF = pd.read_csv(input_file, sep="\t", index_col=0)

# Generate PCA ohne 004_* samples
# Remove all columns containing the string 004_
coresPeaksDF = originalPeaksDF.loc[:,~originalPeaksDF.columns.str.contains('004_', case=False)] 

# Get the total sum for each peak
librarySizeDF = coresPeaksDF.sum(axis=0)

# Principal component analysis
features = coresPeaksDF.columns.tolist()
x        = coresPeaksDF.T
# Standardize the Data: Since PCA yields a feature subspace that maximizes the variance along the axes, 
# it makes sense to standardize the data, especially, if it was measured on different scales. 
# x = StandardScaler().fit_transform(x)

# PCA Projection to 2D
pca   = PCA(n_components=5)
pcs   = pca.fit_transform(x)
pcaDF = pd.DataFrame(data = pcs, columns = ['PC1', 'PC2','PC3', 'PC4','PC5'])


# Add the cellLines information
pcaDF = pd.concat([pcaDF, pd.DataFrame(data=features, columns=['CellLines'])], axis = 1)

# Add groups information
# https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
pcaDF['groups'] = pcaDF["CellLines"].str.split("_", n = 1, expand = True)[0]                                          # 5320, 6075, ...
pcaDF['organs'] = pcaDF["CellLines"].str.split("_", n = 2, expand = True)[1].str.split("-", n = 1, expand = True)[0]  # livmet, lungmet, ppt, ...
pcaDF['groups'] = pcaDF['groups'].apply(np.int)
mapdi ={5320:'m5320', 6075:'m6075', 53631:'m53631', 53646:'m53646'}
pcaDF.replace({"groups": mapdi}, inplace=True)
pcaDF.set_index('CellLines', inplace=True)
pcaDF = pd.concat([pcaDF,librarySizeDF], axis=1)
pcaDF.rename(columns={0:'LibrarySize'}, inplace=True)
pcaDF['CellLines'] = pcaDF.index
pcaDF['sno'] = np.arange(len(pcaDF)) + 1
n = pcaDF['sno'].tolist()

# Perform a Scree Plot of the Principal Components
# A scree plot is like a bar chart showing the size of each of the principal components. 
# It helps us to visualize the percentage of variation captured by each of the principal components
percent_variance = np.round(pca.explained_variance_ratio_* 100, decimals =2)
columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']
screeDF = pd.DataFrame(percent_variance, index=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
screeDF.reset_index(inplace=True)
screeDF.rename(columns={ screeDF.columns[1]: "variance" }, inplace = True)
sns_t = sns.barplot(x='variance',y='index', data=screeDF, palette="Blues_d", orient='h')
show_values_on_bars(sns_t, "h", 0.3)
plt.xlabel('Percentate of Variance Explained')
plt.ylabel('Principal Components')
plt.title('PCA Scree Plot')
# Hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)
# or
ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
# plt.show()
screePlotPdf = "{0}_PCA_Scree_plot.pdf".format(get_file_info(input_file)[3])
screePlotPdf = "{0}_PCA_Scree_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(screePlotPdf,bbox_inches = 'tight')
plt.close('all')

# Visualize 2D Projection
fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='CellLines', style='groups', size='CellLines')
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}_PCA_plot.pdf".format(get_file_info(input_file)[3])
pcaPlotPdf = "{0}_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='organs', s=200, size='LibrarySize')
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}_organs_PCA_plot.pdf".format(get_file_info(input_file)[3])
pcaPlotPdf = "{0}_organs_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='groups')
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}_groups_PCA_plot.pdf".format(get_file_info(input_file)[3])
pcaPlotPdf = "{0}_groups_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')


# Generate the heatmap of top 1% varied peaks
sns.set(font_scale=0.5)
h = sns.clustermap(coresPeaksDF,z_score=0,cmap=sns.diverging_palette(220, 20, n=7),figsize=(10, 20)); 
h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
# plt.show()
heatmapPlotPdf = "{0}_top1pc_heatmap.pdf".format(get_file_info(input_file)[3])
heatmapPlotPdf = "{0}_top1pc_heatmap_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Generate the sample correlation heatmap from top 1% varied peaks
# sns.set(font_scale=0.5)
h = sns.clustermap(coresPeaksDF.corr(),figsize=(10, 10),cmap='gist_heat_r', vmax=1.1, vmin=-0.1)
heatmapPlotPdf = "{0}_sampleCorrelation_top1pc_heatmap.pdf".format(get_file_info(input_file)[3])
heatmapPlotPdf = "{0}_sampleCorrelation_top1pc_heatmap_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

h = sns.clustermap(coresPeaksDF.corr(),figsize=(10, 10),cmap='gist_heat_r', vmax=1.1, vmin=-0.1, annot=True, annot_kws={"size": 5})
heatmapPlotPdf = "{0}_sampleCorrelation_top1pc_heatmap_annotated.pdf".format(get_file_info(input_file)[3])
heatmapPlotPdf = "{0}_sampleCorrelation_top1pc_heatmap_annotated_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')


fig = plt.figure(figsize=(20,15))
ax = fig.add_subplot(1,1,1)
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='groups', s=200)
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
ax=plt.gca()
label_point(pcaDF.PC1, pcaDF.PC2, pcaDF.sno, ax) 
pcaPlotPdf = "{0}_annotated_groups_PCA_plot.pdf".format(get_file_info(input_file)[3])
pcaPlotPdf = "{0}_annotated_groups_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')

######################### USER DEFINED FUNCTIONS ########################
def show_values_on_bars(axs, h_v="v", space=0.4):
    '''
    Two parameters explained: h_v - Whether the barplot is horizontal or vertical. 
                              "h" represents the horizontal barplot, 
                              "v" represents the vertical barplot.
    space - The space between value text and the top edge of the bar. Only works for horizontal mode.
    Source: https://stackoverflow.com/questions/43214978/seaborn-barplot-displaying-values'''
    def _show_on_single_plot(ax):
        if h_v == "v":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height()
                value = p.get_height()
                ax.text(_x, _y, value, ha="center") 
        elif h_v == "h":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() + float(space)
                _y = p.get_y() + p.get_height()
                value = p.get_width()
                ax.text(_x, _y, value, ha="left")

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)

# Generate PCA with 
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.05, point['y']+0.5, str(int(point['val'])))

# CLtr +D +D
#****************************************************************************************************
# Sort input file and remove duplicate lines
f='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_vstCounts.matrix'
(head -n 1 ${f} && tail -n +2 ${f} | sort -u) > tmp && mv tmp ${f}

# Draw heatmap with clusterIDs and save the clusterIDs in a separate file
R
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(matrixStats))

## Get the input data
cat("- Reading input file ...\n")
inputfile              <- '/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/cream/merged/5320_53631_53646_6075_merge_master_peaks_COREs_vstCounts.matrix'
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
numClusters          <- 10
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
finHeatmap      <- pheatmap(scaled_df, filename=pdffile, cutree_row=numClusters, show_rownames=F, annotation_row=annClust, height=30, width=10)



# bam to bigwig using deeptools bamCoverage
for f in /media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/bams/trimmed/*.bam
do
 echo ${f}
 bname=$(basename ${f} .bam)
 od=$(dirname ${f})/bigwig
 bamCoverage -b ${f} -o ${od}/${bname}.bw -p 64
 echo ""
done

