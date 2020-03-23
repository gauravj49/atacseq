# PWD
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# 1) Get the peaks in the cytokine region
peakMatrix="/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_vstCounts.tab"

# 2.2) Get bed files for each sample
od="$(dirname ${peakMatrix})/samplesBed"
mkdir -p ${od}
header=$(head -1 ${peakMatrix})
numCols=$(awk -F'\t' '{print NF; exit}' ${peakMatrix})
for i in {4..${numCols}};
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

# 2.3) Get the mean score of atacseq peaks for the cytokine region
pd="$(dirname ${peakMatrix})/cytokinePeaks"; mkdir -p ${pd}
ib="/home/rad/users/gaurav/projects/seqAnalysis/atacseq/docs/Interleukine_list.bed"
sort -k1,1 -k2n ${ib} -o ${ib}
for bbed in ${od}/*.bed;
do
 echo ${bbed}
 bedtools map -a ${ib} -b ${bbed} -c 4 -o mean > ${pd}/$(basename ${bbed} .bed)_cytokine_peaks_mean.bed
 # Add header to the file
 header="chr\tstart\tend\tname\t$(basename ${bbed} _atacseq_mm_se_rmdup_peaks.bed)_mean"
 sed -i "1i${header}" ${pd}/$(basename ${bbed} .bed)_cytokine_peaks_mean.bed
done

# Merge the individual mean peaks into a dataframe
ipython
from scipy.stats import zscore

# Read and merge into df
pkDF = pd.concat([pd.read_csv(f, sep='\t').set_index(['chr', 'start', 'end', 'name']) for f in glob.glob('/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/cytokinePeaks/*.bed')],axis=1).reset_index()
pkDF.replace('.',0, inplace=True)

# Save to 
output_file = '/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_Interleukine_merge_peaks_averageReads.txt'
pkDF.to_csv(output_file, sep='\t', index = False, float_format='%.2g')

# Plot the heatmap
peaksDF = pkDF.copy()
peaksDF.drop(columns = ['chr','start','end'], inplace = True)
peaksDF.set_index('name',inplace=True)
peaksDF = peaksDF.loc[~(peaksDF==0).all(axis=1)] # Remove rows with all 0s
peaksDF = peaksDF.astype('float') # Conver all columns to float
# Remove _peaks.bed_mean from column names
peaksDF = peaksDF.rename(columns={col: col.split('_peaks')[0] for col in peaksDF.columns})

# sns.clustermap(peaksDF, z_score=1, cmap='RdBu_r')
# heatmapPlotPdf = "{0}_heatmap.pdf".format(get_file_info(output_file)[3])
# plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
# plt.close('all')

# Rearrange peaksDF columns
new_columns_order = ['SB02_DP_CAGAGAGG','SB01_DP_TAAGGCGA','GS250_DP_CTCTCT','GS251_DP_CAGAGA','GS244_12h_NKT_TAAGGC','GS245_12h_NKT_CGTACT','GS394_24h_NKT_TAAGGCGA','SB04_24h_NKT_CGAGGCTG','GS395_24h_NKT_TAAGGCGA','SB03_24h_NKT_CGAGGCTG','SB06_36h_NKT_AAGAGGCA','GS396_36h_NKT_CGTACTAG','GS397_36h_NKT_CGTACTAG','GS398_48h_NKT_AGGCAGAA','SB09_48h_NKT_CGTACTAG','SB10_48h_NKT_CGTACTAG','SB13_60h_NKT_AGGCAGAA','GS399_60h_NKT_TCCTGAGC','SB12_60h_NKT_CGTACTAG','GS400_72h_NKT_GGACTCCT','SB17_72h_NKT_GTCGTGAT','SB18_72h_NKT_CTCTCTAC','SB16_72h_NKT_CTCTCTAC','SB15_72h_NKT_TCCTGAGC','GS401_84h_NKT_TCCTGAGC','SB21_84h_NKT_GGACTCCT','SB22_96h_NKT_TAGGCATG','GS402_96h_NKT_TAGGCATG','GS247_5d_NKT_TCCTGA','GS246_5d_NKT_AGGCAG','GS403_NKT1_CTCTCTAC','GS248_NKT1_GGACTC','GS249_NKT1_TAGGCA']

# NOTE: here the sample (GS400_71h_*) belongs to 72h group so using that in the groups annotation below
groups = ['DP','DP','DP','DP','12h','12h','24h','24h','24h','24h','36h','36h','36h','48h','48h','48h','60h','60h','60h','72h','72h','72h','72h','72h','84h','84h','96h','96h','120h','120h','NKT1','NKT1','NKT1']

# Rearrange columns to new order
peaksDF = peaksDF[new_columns_order]
peaksDF = peaksDF.astype('float') # Conver all columns to floa

# Plot the dataframe as heatmap
sns.clustermap(peaksDF, z_score=1, cmap='RdBu_r', col_cluster=False)
heatmapPlotPdf = "{0}_heatmap_rearranged.pdf".format(get_file_info(output_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Add groups to the peaksDF
peaksDF.loc['groups'] = groups

# Get the meanPeaksDF for each group
meanPeaksDF = pd.DataFrame()
peaksGroups = peaksDF.T.groupby('groups')

for groupName, groupDF in peaksGroups:
        groupDF.drop('groups', inplace=True, axis=1)
        meanPeaksDF = pd.concat([meanPeaksDF, pd.DataFrame(groupDF.mean(axis=0), columns=['mean_{0}'.format(groupName)])], axis=1)

# Rename the index
meanPeaksDF.index.name = 'featureid'

# Rearrange columns to new order
mean_col_order = ['mean_DP','mean_12h','mean_24h','mean_36h','mean_48h','mean_60h','mean_72h','mean_84h','mean_96h','mean_120h','mean_NKT1']
meanPeaksDF    = meanPeaksDF[mean_col_order]
meanPeaksDF     = meanPeaksDF.astype('float')

# Calculate the zscore of the dataframe
# The peaksDF.apply(zscore, axis=1) returns an array 
# https://stackoverflow.com/questions/35491274/pandas-split-column-of-lists-into-multiple-columns
from scipy.stats import zscore
peaksDF.drop('groups', inplace=True)
peaksDF = peaksDF.astype('float')
zscorePeaksDF     = pd.DataFrame(pd.DataFrame(peaksDF.apply(zscore, axis=1))[0].values.tolist(), columns=peaksDF.columns.tolist(), index=peaksDF.index)
zscoreMeanPeaksDF = pd.DataFrame(pd.DataFrame(meanPeaksDF.apply(zscore, axis=1))[0].values.tolist(), columns=meanPeaksDF.columns.tolist(), index=meanPeaksDF.index)

# Plot the heatmap for all replicates separately for all time points
g = sns.clustermap(zscorePeaksDF, cmap='RdBu_r', col_cluster=False, figsize=(10,8)); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0); cax = plt.gcf().axes[-1]; cax.tick_params(labelsize=5);
heatmapPlotPdf = "{0}_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')

# Plot the heatmap for mean of replicates for all time points
g = sns.clustermap(zscoreMeanPeaksDF, cmap='RdBu_r', col_cluster=False, figsize=(10,8)); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0); cax = plt.gcf().axes[-1]; cax.tick_params(labelsize=5);
heatmapPlotPdf = "{0}_of_replicates_for_timepoints_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')

# Save the meanPeaksDf
meanPeaksDF['meanFeatureid'] = meanPeaksDF.mean(axis=1)
meanPeaksDF.sort_values(inplace=True, by='meanFeatureid')
output_file = '/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_Interleukine_merge_peaks_mean_of_replicates_for_timepoints.txt'
meanPeaksDF.to_csv(output_file, sep='\t', index = True, float_format='%.2g')

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

####################### USER DEFINIED FUNCTIONS ###########################
def fixedWidthClusterMap(dataFrame, cellSizePixels=75):
        '''
        Make multiple plots where the cell sizes are exactly identical
        '''

        # Calulate the figure size, this gets us close, but not quite to the right place
        dpi = matplotlib.rcParams['figure.dpi']
        marginWidth = matplotlib.rcParams['figure.subplot.right']-matplotlib.rcParams['figure.subplot.left']
        marginHeight = matplotlib.rcParams['figure.subplot.top']-matplotlib.rcParams['figure.subplot.bottom']
        Ny,Nx = dataFrame.shape
        figWidth = (Nx*cellSizePixels/dpi)/0.8/marginWidth
        figHeigh = (Ny*cellSizePixels/dpi)/0.8/marginHeight

        # # do the actual plot
        # grid = sns.clustermap(dataFrame, figsize=(figWidth, figHeigh))
        grid = sns.clustermap(dataFrame, z_score=1, cmap='RdBu_r', col_cluster=False, vmin=-1.7, vmax=1.7);

        # calculate the size of the heatmap axes
        axWidth = (Nx*cellSizePixels)/(figWidth*dpi)
        axHeight = (Ny*cellSizePixels)/(figHeigh*dpi)

        # resize heatmap
        ax_heatmap_orig_pos = grid.ax_heatmap.get_position()
        grid.ax_heatmap.set_position([ax_heatmap_orig_pos.x0, ax_heatmap_orig_pos.y0, axWidth, axHeight])

        # resize dendrograms to match
        ax_row_orig_pos = grid.ax_row_dendrogram.get_position()
        grid.ax_row_dendrogram.set_position([ax_row_orig_pos.x0, ax_row_orig_pos.y0, ax_row_orig_pos.width, axHeight])
        ax_col_orig_pos = grid.ax_col_dendrogram.get_position()
        grid.ax_col_dendrogram.set_position([ax_col_orig_pos.x0, ax_heatmap_orig_pos.y0+axHeight, axWidth, ax_col_orig_pos.height])
        return grid # return ClusterGrid object

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

