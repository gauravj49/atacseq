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
#****************************************************************************************************
pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

input_file  = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.tab"
output_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.matrix"
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

# Create binary peaks flag annotation matrix
peaksTabFile='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.tab'
peaksAnnFile='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_annotation.txt'
peakFilesDir='/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/peakFiles'
tempPkAnnDir="${peakFilesDir}/tempPkAnn"; mkdir -p ${tempPkAnnDir}
colstofilter="," # Column numbers to be extracted at the end
i=1              #
header="Chrom\tStart\tEnd\tPeakID\t"        
for p in ${peakFilesDir}/*.narrowPeak;
do
 bname=$(basename ${p} _R1_001_rmdup_peaks.narrowPeak)
 outfile=${tempPkAnnDir}/${bname}.bed
 intersectBed -a ${peaksTabFile} -b ${p} -c | awk '{print $1,$2,$3,$4,$NF}' > ${outfile}
 colstofilter=$(echo "${colstofilter}$((i*5)),");
 header=$(echo -e "${header}${bname}\t")
 echo "${i}) ${bname}: ${colstofilter}"
 echo ""
 i=$((i+1));
done

# Remove last "," from the string
colstofilter=$(echo ${colstofilter}|sed 's/,$//')
# ,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90

# Paste output the file with multiple delimiters
# Columns from one files are space separated and multiple files are separated by space
# Using sed to convert the tabs to spaces and then using cut to get the final columns
paste *.bed| sed 's/\t/ /g' | cut -d ' ' --output-delimiter=$'\t' -f1-4${colstofilter}> ${peaksAnnFile}

# Add the header to the file
sed  -i "1i${header}" ${peaksAnnFile}

############# DOCS #############
# From the pca, it seems
# PC1 = Morphology i.e. (c2a, c2b and c2c)
# PC2 = Library Size which is suppose to be organs
# PC3 = Mouse
################################
# Get the top 1% most varying genes and plot the pca
ipython
#****************************************************************************************************
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

# Input and output files
input_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.matrix"
rankoutput_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/5320_53631_53646_6075_merge_master_peaks_rawCounts.ranks"

# Import data into the dataframe
peaksDF    = pd.read_csv(input_file, sep="\t", index_col=0)

# Get the total sum for each peak
librarySizeDF = peaksDF.sum(axis=0)

# Get the variance for each peak
peaksDF['variance'] = peaksDF.var(axis=1)

# Rank the variance (expressed as percentile rank)
peaksDF['pctRank']  = peaksDF['variance'].rank(pct=True)

# Save the variance and ranks genes
peaksDF[['variance','pctRank']].to_csv(rankoutput_file, index=True, header=True, sep="\t", float_format='%.2f')

# Get the filtered dataframe by taking top 1% of ranked peaks
topPeaksDF = peaksDF[peaksDF['pctRank'] >= 0.99]

# Get the dimensions of original and top peaks df
peaksDF.shape    # peaksDF.shape
topPeaksDF.shape # (1600, 20)

# Drop the variance pctRnak columns from the top ranked genes
topPeaksDF.drop(columns=['variance','pctRank'], inplace=True)

# Principal component analysis
features = topPeaksDF.columns.tolist()
x        = topPeaksDF.T
# Standardize the Data: Since PCA yields a feature subspace that maximizes the variance along the axes, 
# it makes sense to standardize the data, especially, if it was measured on different scales. 
x = StandardScaler().fit_transform(x)

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
# screePlotPdf = "{0}_PCA_Scree_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
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
# pcaPlotPdf = "{0}_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')

pcaDF.set_index('CellLines', inplace=True)
pcaDF = pd.concat([pcaDF,librarySizeDF], axis=1)
pcaDF.rename(columns={0:'LibrarySize'}, inplace=True)
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
# pcaPlotPdf = "{0}_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
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
# pcaPlotPdf = "{0}_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')


# Generate the heatmap of top 1% varied peaks
sns.set(font_scale=0.5)
h = sns.clustermap(topPeaksDF,z_score=0,cmap=sns.diverging_palette(220, 20, n=7),figsize=(10, 20)); 
h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
# plt.show()
heatmapPlotPdf = "{0}_top1pc_heatmap.pdf".format(get_file_info(input_file)[3])
heatmapPlotPdf = "{0}_top1pc_heatmap_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Generate the sample correlation heatmap from top 1% varied peaks
# sns.set(font_scale=0.5)
h = sns.clustermap(topPeaksDF.corr(),figsize=(10, 10),cmap='gist_heat_r', vmax=1.1, vmin=-0.1)
# plt.show()
heatmapPlotPdf = "{0}_sampleCorrelation_top1pc_heatmap.pdf".format(get_file_info(input_file)[3])
heatmapPlotPdf = "{0}_sampleCorrelation_top1pc_heatmap_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')




# Generate PCA ohne 004_* samples
# Remove all columns containing the string 004_
originalPeaksDF=peaksDF.copy()
peaksDF = peaksDF.loc[:,~peaksDF.columns.str.contains('004_', case=False)] 
# peaksDF.shape =  (159990, 17)

# Follow the same code from above to generate the ohneOutliers plots
# Make sure to change the file names







# Pairwise similarity matrices analysis
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
#****************************************************************************************************
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
#****************************************************************************************************

# Plot similarity heatmap
ipython
#****************************************************************************************************
sns.set(font_scale=0.3)
peaksDF = pd.read_csv("/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_distance.matrix", sep="\t", index_col=0) 
g = sns.clustermap(peaksDF, vmax=0.02, cmap="ocean"); 
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 5)
plt.savefig("/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/jaccard_pairwise_peaks_comparisons_heatmap.pdf")               
#****************************************************************************************************

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
