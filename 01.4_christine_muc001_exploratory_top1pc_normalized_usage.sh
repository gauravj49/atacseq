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

pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

# Input and output files
input_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/5320_53631_53646_6075_merge_master_peaks_vstCounts.matrix"
rankoutput_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/5320_53631_53646_6075_merge_master_peaks_vstCounts.ranks"

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
# x = StandardScaler().fit_transform(x)

# PCA Projection to 2D
pca   = PCA(n_components=5)
pcs   = pca.fit_transform(x)
pcs   = pca.fit(x)
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
