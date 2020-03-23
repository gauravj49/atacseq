# PWD
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

############# DOCS #############
################################
# Get the top 1% most varying genes and plot the pca
ipython
#****************************************************************************************************
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

# Input and output files
input_file      = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_vstCounts.matrix"
rankoutput_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_vstCounts.ranks"
read_cutoff     = 5
sample_pc       = 75

# Import data into the dataframe
peaksDF    = pd.read_csv(input_file, sep="\t", index_col=0)

# Get the total sum for each peak
librarySizeDF = peaksDF.sum(axis=0)

# Filter peaks sum of 20 in all the samples
origPeaksDF = peaksDF.copy()

# Filter peaks with sum of peaks less than 75 for all samples
peaksDF = peaksDF[peaksDF.sum(axis=1) >= 75]

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

# Save the topPeaks files in tab separated files
# tabPeaksDF = topPeaksDF.copy()
tabPeaksDF = peaksDF.copy()

# Get the chr, start and end columns
tabPeaksDF['tabDFSlice'] = tabPeaksDF.index.values.tolist()
tabDFSlice               = tabPeaksDF['tabDFSlice'].str.split("_", n = 3, expand = True)
tabPeaksDF['chr']        = tabDFSlice[0]
tabPeaksDF['start']      = tabDFSlice[1]
tabPeaksDF['end']        = tabDFSlice[2]

# Rearrange columns to new order
tabPeaksDF.drop(columns  =['variance','pctRank'], inplace=True)
tabPeaksDF.rename(columns={'tabDFSlice':'PeakID', 'chr':'PeakChr', 'start':'PeakStart', 'end':'PeakEnd'}, inplace=True)
new_columns_order        = ['PeakChr', 'PeakStart', 'PeakEnd','PeakID', 'TransPB-Bcell-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-Bcell-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-CD4-SinglePositive-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-CD4-SinglePositive-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-CD8-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-CD8-SinglePositive-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-CLP-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-CLP-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-DN2a-Tcell-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-DN2a-Tcell-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-DN2b-Tcell-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-DN2b-Tcell-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-DN3-Tcell-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-DN3-Tcell-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-DN4-Tcell-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-DN4-Tcell-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-DP-Tcell-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-DP-Tcell-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-ETP-Tcell-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-ETP-Tcell-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-HSC-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-HSC-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-MPP-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-MPP-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup', 'TransPB-NK-BoneMarrow-Rep1_mm_atacseq_pe_R1_rmdup', 'TransPB-NK-BoneMarrow-Rep2_mm_atacseq_pe_R1_rmdup',]
tabPeaksDF               = tabPeaksDF[new_columns_order]

# Save the DF as tab separated file
taboutput_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_vstCounts.tab"
tabPeaksDF.to_csv(taboutput_file, sep='\t', index = False, float_format='%.2g')

# Save the top1pcDF as tab separated file
tabPeaksDF.set_index('PeakID', inplace=True)
tabPeaksDF['PeakID'] = tabPeaksDF.index
bedPeaksDF = tabPeaksDF.iloc[:,:4].copy()
#default inner join
top1pcbedDF = pd.merge(bedPeaksDF, topPeaksDF, left_index=True, right_index=True)

top1pctaboutput_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/cis/tALLcellLineMm_cis_all_merge_master_peaks_vstCounts_top1pc.tab"
tabPeaksDF.to_csv(top1pctaboutput_file, sep='\t', index = False, float_format='%.2g')



# Principal component analysis
features = topPeaksDF.columns.tolist()
x        = topPeaksDF.T
# Standardize the Data: Since PCA yields a feature subspace that maximizes the variance along the axes, 
# it makes sense to standardize the data, especially, if it was measured on different scales. 
# x = StandardScaler().fit_transform(x)

# PCA Projection to 2D
pca   = PCA(n_components=5)
pcs   = pca.fit_transform(x)
pcaDF = pd.DataFrame(data = pcs, columns = ['PC1', 'PC2','PC3', 'PC4','PC5'])

# Add the cellLines information
pcaDF = pd.concat([pcaDF, pd.DataFrame(data=features, columns=['CellLines'])], axis = 1)

# # Add batch information
# pcaDF['batch'] = sum([['B1'] * 8,['B2'] * 10,['B3'] * 15], [])
# pcaDF = pd.concat([pcaDF, pd.DataFrame(data=timeline, columns=['timeline'])], axis = 1)

# Add groups information
# https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
pcaDF['timeline'] = pcaDF["CellLines"].str.split("_", n = 2, expand = True)[0].str.replace('TransPB-', '') # Bcell-BoneMarrow-Rep1
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
plt.savefig(screePlotPdf,bbox_inches = 'tight')
plt.close('all')

# Visualize 2D Projection
filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
# g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='CellLines', style='timeline', size='CellLines', markers=filled_markers)
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='CellLines', size='CellLines', markers=filled_markers)
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}_PCA_plot.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
# g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='timeline', s=200, size='LibrarySize')
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='timeline', s=200)
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}_timeline_PCA_plot.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')



# UMAP
import umap
reducer   = umap.UMAP(n_components=5)
embedding = reducer.fit_transform(topPeaksDF.T)
umapDF = pd.DataFrame(data = embedding, columns = ['UMap1', 'UMap2','UMap3', 'UMap4','UMap5'])
timeline = [x.split('_')[0].replace('TransPB-', '') for x in features]
umapDF = pd.concat([umapDF, pd.DataFrame(data=timeline, columns=['timeline'])], axis = 1)
g = sns.scatterplot(x='UMap1', y='UMap2', data=umapDF, hue='timeline', s=100, alpha=0.8, palette='viridis')
plt.gca().set_aspect('equal', 'datalim')
plt.gca().set_title('UMAP projection of Top 1% variance ranked peaks', fontsize = 20)
umapPlotPdf = "{0}_timeline_UMAP_plot.pdf".format(get_file_info(input_file)[3])
plt.savefig(umapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Generate the heatmap of top 1% varied peaks
sns.set(font_scale=0.5)
h = sns.clustermap(topPeaksDF,z_score=0,cmap=sns.diverging_palette(220, 20, n=7),figsize=(10, 20)); 
h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
# plt.show()
heatmapPlotPdf = "{0}_top1pc_heatmap.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Generate the sample correlation heatmap from top 1% varied peaks
# sns.set(font_scale=0.5)
h = sns.clustermap(topPeaksDF.corr(),figsize=(10, 10),cmap='gist_heat_r', vmax=1.1, vmin=-0.1)
# plt.show()
heatmapPlotPdf = "{0}_sampleCorrelation_top1pc_heatmap.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Generate PCA with 
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.05, point['y']+0.5, str(int(point['val'])))
        
fig = plt.figure(figsize=(20,15))
ax = fig.add_subplot(1,1,1)
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='timeline', s=200)
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
ax=plt.gca()
label_point(pcaDF.PC1, pcaDF.PC2, pcaDF.sno, ax) 
pcaPlotPdf = "{0}_annotated_timeline_PCA_plot.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')

# Save the cellLines information in a file
pcaCellLinesFile = "{0}_annotated_groups_PCA.txt".format(get_file_info(input_file)[3])
pcaDF[['sno','CellLines']].to_csv(pcaCellLinesFile, sep="\t", index=None)

####################### USER DEFINIED FUNCTIONS ###########################
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




