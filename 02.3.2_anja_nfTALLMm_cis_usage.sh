# PWD
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

############# DOCS #############
################################
ipython
#****************************************************************************************************
import datatable as dt
from datatable import *
import scanpy as sc
import time
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.lines import Line2D
filled_markers = itertools.cycle(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
fillstyles     = ('full', 'left', 'right', 'bottom', 'top', 'none')


pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

####################### USER DEFINIED FUNCTIONS ###########################
# Generate PCA with 
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.05, point['y']+0.5, str(int(point['val'])))

def show_values_on_bars(axs, h_v="v", space=0.4):
        '''https://stackoverflow.com/questions/43214978/seaborn-barplot-displaying-values'''
        def _show_on_single_plot(ax):
                if h_v == "v":
                        for p in ax.patches:
                                _x = p.get_x() + p.get_width() / 2
                                _y = p.get_y() + p.get_height()
                                value = int(p.get_height())
                                ax.text(_x, _y, value, ha="center") 
                elif h_v == "h":
                        for p in ax.patches:
                                _x = p.get_x() + p.get_width() + float(space)
                                _y = p.get_y() + p.get_height()
                                value = int(p.get_width())
                                ax.text(_x, _y, value, ha="left")

                if isinstance(axs, np.ndarray):
                        for idx, ax in np.ndenumerate(axs):
                                _show_on_single_plot(ax)
                else:
                        _show_on_single_plot(axs)

################################################################################
# Input and output files
input_file = "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_analysis_consensus_peaks_vstCounts.txt"
annot_file = "/media/rad/HDD1/atacseq/anja/nfTALLMm/analysis/nfTALLMm_sample_annotation.txt"
projName   = "nfTALLMm"
bname      = "{0}_mLcP_vst".format(projName) # mLcP = merged_library_consensus_peaks
outputDir  = "{0}/cis".format(get_file_info(input_file)[0]); create_dir(outputDir)
plotsDir   = "{0}/plots".format(outputDir); create_dir(plotsDir)
dataDir    = "{0}/data".format(outputDir); create_dir(dataDir)
sample_pc  = 75

# Import data into the datatable
start_time = time.time()
peaksDF    = pd.read_csv(input_file, sep="\t", index_col=3)
annotDF    = pd.read_csv(annot_file, sep="\t", index_col=0)
print("%s seconds" % (time.time() - start_time))
# 0.0713660717010498 seconds

# Rearrange columns to new order
new_columns_order        = ['PeakChrom', 'PeakStart', 'PeakEnd','HSCBoneMarrow_R1','HSCBoneMarrow_R2','MPPBoneMarrow_R1','MPPBoneMarrow_R2','CLPBoneMarrow_R1','CLPBoneMarrow_R2','ETPTcell_R1','ETPTcell_R2','DN2aTcell_R1','DN2aTcell_R2','DN2bTcell_R1','DN2bTcell_R2','DN3Tcell_R1','DN3Tcell_R2','DN4Tcell_R1','DN4Tcell_R2','DPTcell_R1','DPTcell_R2','CD4SinglePositive_R1','CD4SinglePositive_R2','CD8SinglePositive_R1','CD8SinglePositive_R2','BcellBoneMarrow_R1','BcellBoneMarrow_R2','NKBoneMarrow_R1','NKBoneMarrow_R2']
peaksDF = peaksDF[new_columns_order]

# Drop the genomic ranges columns from the peaksDT
peaksDF.drop(columns = ['PeakChrom','PeakStart','PeakEnd'], inplace=True)

# Get the total sum for each peak
librarySizeDF = peaksDF.sum(axis=0)

# Filter peaks sum of 20 in all the samples
origPeaksDF = peaksDF.copy()

# Filter peaks with sum of peaks less than 75 for all samples
peaksDF = peaksDF[peaksDF.sum(axis=1) >= sample_pc]

# Get the variance for each peak
peaksDF['variance'] = peaksDF.var(axis=1)

# Rank the variance (expressed as percentile rank)
peaksDF['pctRank']  = peaksDF['variance'].rank(pct=True)

# Save the variance and ranks genes
peaksDF[['variance','pctRank']].to_csv("{0}/{1}.ranks".format(dataDir,bname), index=True, index_label = 'PeakID', header=True, sep="\t", float_format='%.4g')

# Get the filtered dataframe by taking top 5% of ranked peaks
hvPeaksDF = peaksDF[peaksDF['pctRank'] >= 0.95]

# Get the dimensions of original and top peaks df
peaksDF.shape    # peaksDF.shape
hvPeaksDF.shape # (2259, 28)

# Drop the variance pctRnak columns from the top ranked genes
peaksDF.drop  (columns=['variance','pctRank'], inplace=True)
hvPeaksDF.drop(columns=['variance','pctRank'], inplace=True)

# Save the highly variable peaks to a file
hvPeaksDF.to_csv("{0}/{1}_hvPeaks5pc.txt".format(dataDir,bname), index=True, index_label = 'PeakID', header=True, sep="\t", float_format='%.4g')

######################## PLOTS ###########################################
# 2) Principal component analysis
# Standardize the Data: Since PCA yields a feature subspace that maximizes the variance along the axes, 
# it makes sense to standardize the data, especially, if it was measured on different scales. 
# x = StandardScaler().fit_transform(x)
x       = peaksDF.T
pca     = PCA(n_components=5)
pcs     = pca.fit_transform(x)
pcaDF   = pd.DataFrame(data = pcs, columns = ['PC1', 'PC2','PC3', 'PC4','PC5'])

# Add index to pcaDF
pcaDF['pcaIndex'] = peaksDF.columns.tolist()

# Add groups information
# https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
pcaDF.set_index('pcaIndex', inplace=True)
pcaDF = pd.concat([pcaDF,librarySizeDF], axis=1)
pcaDF.rename(columns={0:'LibrarySize'}, inplace=True)
pcaDF = pd.concat([pcaDF,annotDF], axis=1)
pcaDF['sno'] = np.arange(len(pcaDF)) + 1
n = pcaDF['sno'].tolist()

# 2.1) Perform a Scree Plot of the Principal Components
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
screePlotPdf = "{0}/01_{1}_PCA_Scree_plot.pdf".format(plotsDir, bname)
plt.savefig(screePlotPdf,bbox_inches = 'tight')
plt.close('all')

# 2.2) Visualize 2D Projection
filled_markers = itertools.cycle(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
pHue = 'CellLine'; pSize = 'Study'; pStyle = 'CellLine'
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue=pHue, size=pSize, style=pStyle, markers=filled_markers)
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('VST normalized peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}/02_{1}_PCA_col{2}_size{3}_style{4}.pdf".format(plotsDir, bname, pHue, pSize, pStyle)
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
plt.close('all')

# 2.3) PCA with annoatation
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
pHue = 'CellLine'; pSize = 'Study'; pStyle = 'CellLine'
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue=pHue, size=pSize, style=pStyle, markers=filled_markers, linewidth=0)
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 5% variance ranked peaks PCA', fontsize = 20)
ax=plt.gca()
label_point(pcaDF.PC1, pcaDF.PC2, pcaDF.sno, ax) 
pcaPlotPdf = "{0}/03_{1}_annotated_PCA_col{2}_size{3}_style{4}.pdf".format(plotsDir, bname, pHue, pSize, pStyle)
plt.savefig(pcaPlotPdf, bbox_inches = 'tight')
plt.close('all')
# Save the CellLine information in a file
pcaCellLineFile = "{0}/03_{1}_annotated_groups_PCA_col{2}_size200_style{4}.txt".format(plotsDir, bname, pHue, pSize, pStyle)
# pcaDF[['sno','SampleName']].to_csv(pcaCellLineFile, sep="\t", index=None)
pcaDF[['sno']].to_csv(pcaCellLineFile, sep="\t", index=True, index_label='SampleName')

# 2.4) UMAP
import umap
reducer   = umap.UMAP(n_components=5, random_state=2105, n_neighbors=8)
embedding = reducer.fit_transform(peaksDF.T)
umapDF = pd.DataFrame(data = embedding, columns = ['UMAP1', 'UMAP2','UMAP3', 'UMAP4','UMAP5'])
# Add groups information
umapDF['umapIndex'] = peaksDF.columns.tolist()
umapDF.set_index('umapIndex', inplace=True)
umapDF = pd.concat([umapDF,librarySizeDF], axis=1)
umapDF.rename(columns={0:'LibrarySize'}, inplace=True)
umapDF = pd.concat([umapDF,annotDF], axis=1)
umapDF['sno'] = np.arange(len(umapDF)) + 1
n = umapDF['sno'].tolist()
fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot(1,1,1) 
pHue = 'CellLine'; pSize = 'Study'; pStyle = 'CellLine'
g = sns.scatterplot(x='UMAP1', y='UMAP2', data=umapDF, hue=pHue, style=pStyle, palette=sns.color_palette("hls",umapDF[pHue].value_counts().shape[0]), s=100, edgecolor='k', linewidth=0.05, alpha=0.9, markers=filled_markers)
g.legend(loc='lower left', bbox_to_anchor=(1.05, 0.5), ncol=3)
plt.gca().set_aspect('equal', 'datalim')
plt.gca().set_title('UMAP projection of VST normalized peaks', fontsize = 20)
umapPlotPdf = "{0}/04_{1}_UMAP_col{2}_size200_style{3}.pdf".format(plotsDir, bname, pHue, pStyle)
plt.tight_layout(); plt.savefig(umapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# 2.5)Generate the sample correlation heatmap from top 5% varied peaks
h = sns.clustermap(peaksDF.corr(),figsize=(10, 10),cmap='gist_heat_r')
# plt.show()
heatmapPlotPdf = "{0}/05_{1}_sampleCorrelation_heatmap.pdf".format(plotsDir, bname)
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# 2.6) Generate the heatmap of top 5% varied peaks
sns.set(font_scale=0.5)
h = sns.clustermap(hvPeaksDF,z_score=0,cmap=sns.diverging_palette(220, 20, n=7),figsize=(10, 20)); 
h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
# plt.show()
heatmapPlotPdf = "{0}/06_{1}_HighlyVariablePeaks5pc_heatmap.pdf".format(plotsDir, bname)
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')
