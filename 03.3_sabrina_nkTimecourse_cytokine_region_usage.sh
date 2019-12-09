# PWD
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# 1) Get the tab seaprated raw counts of the merged peaks for all the samples
# Using deeptools
multiBamSummary BED-file --BED /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks.bed --bamfiles /media/rad/HDD1/atacseq/sabrina/nkTimecourse/mapping/bams/trimmed/*.bam --smartLabels -out /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.npz --outRawCounts /media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab -p 64

# 2) Get the peaks in the cytokine region
peakMatrix="/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_all_merge_master_peaks_rawCounts.tab"
# # 2.1) Remove single quotes from the header
# sed -i "s/'//g" ${peakMatrix}

# 2.2) Get bed files for each sample
od="$(dirname ${peakMatrix})/samplesBed"
mkdir -p ${od}
header=$(head -1 ${peakMatrix})
for i in {4..36};
do
 bname=$(cut -f${i} <<< ${header})
 echo ${bname}
 cut -f1-3,${i} ${peakMatrix} > ${od}/${bname}_peaks.bed
 # Sort the bed file 
 sort -k1,1 -k2n ${od}/${bname}_peaks.bed -o ${od}/${bname}_peaks.bed
done

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

# Read and merge into df
pkDF = pd.concat([pd.read_csv(f, sep='\t').set_index(['chr', 'start', 'end', 'name']) for f in glob.glob('/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/cytokinePeaks/*.bed')],axis=1).reset_index()

# Save to 
output_file = '/media/rad/HDD1/atacseq/sabrina/nkTimecourse/analysis/nkTimecourse_Interleukine_merge_peaks_mean.txt'
pkDF.to_csv(output_file, sep='\t', index = False, float_format='%.2g')

# Plot the heatmap
peaksDF = pkDF.copy()
peaksDF.drop(columns = ['chr','start','end'], inplace = True)
peaksDF.set_index('name',inplace=True)
sns.clustermap(peaksDF, z_score=1, cmap='RdBu_r')
heatmapPlotPdf = "{0}_heatmap.pdf".format(get_file_info(output_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Rearrange peaksDF columns
new_columns_order = ['SB02_DP_CAGAGAGG_mean','SB01_DP_TAAGGCGA_mean','GS250_DP_CTCTCT_mean','GS251_DP_CAGAGA_mean','GS244_12h_NKT_TAAGGC_mean','GS245_12h_NKT_CGTACT_mean','GS394_24h_NKT_TAAGGCGA_mean','SB04_24h_NKT_CGAGGCTG_mean','GS395_24h_NKT_TAAGGCGA_mean','SB03_24h_NKT_CGAGGCTG_mean','SB06_36h_NKT_AAGAGGCA_mean','GS396_36h_NKT_CGTACTAG_mean','GS397_36h_NKT_CGTACTAG_mean','GS398_48h_NKT_AGGCAGAA_mean','SB09_48h_NKT_CGTACTAG_mean','SB10_48h_NKT_CGTACTAG_mean','SB13_60h_NKT_AGGCAGAA_mean','GS399_60h_NKT_TCCTGAGC_mean','SB12_60h_NKT_CGTACTAG_mean','GS400_71h_NKT_GGACTCCT_mean','SB17_72h_NKT_GTCGTGAT_mean','SB18_72h_NKT_CTCTCTAC_mean','SB16_72h_NKT_CTCTCTAC_mean','SB15_72h_NKT_TCCTGAGC_mean','GS401_84h_NKT_TCCTGAGC_mean','SB21_84h_NKT_GGACTCCT_mean','SB22_96h_NKT_TAGGCATG_mean','GS402_96h_NKT_TAGGCATG_mean','GS247_5d_NKT_TCCTGA_mean','GS246_5d_NKT_AGGCAG_mean','GS403_NKT1_CTCTCTAC_mean','GS248_NKT1_GGACTC_mean','GS249_NKT1_TAGGCA_mean']

# NOTE: here the sample (GS400_71h_*) belongs to 72h group so using that in the groups annotation below
groups = ['DP','DP','DP','DP','12h','12h','24h','24h','24h','24h','36h','36h','36h','48h','48h','48h','60h','60h','60h','72h','72h','72h','72h','72h','84h','84h','96h','96h','120h','120h','NKT1','NKT1','NKT1']

# Rearrange columns to new order
peaksDF = peaksDF[new_columns_order]

# Add groups to the peaksDF
peaksDF.loc['groups'] = groups

# Plot the dataframe as heatmap
sns.clustermap(peaksDF, z_score=1, cmap='RdBu_r', col_cluster=False)
heatmapPlotPdf = "{0}_heatmap_rearranged.pdf".format(get_file_info(output_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Plot heatmap ohne TNF
peaksOhneTNFDF = peaksDF.drop('TNF')
sns.clustermap(peaksOhneTNFDF, z_score=1, cmap='RdBu_r', col_cluster=False)
heatmapPlotPdf = "{0}_heatmap_rearranged_ohne_TNF.pdf".format(get_file_info(output_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Plot heatmap ohne TNF and Tgfb1
peaksOhneTNFTgfb1DF = peaksDF.drop(['TNF', 'Tgfb1'])
sns.clustermap(peaksOhneTNFTgfb1DF, z_score=1, cmap='RdBu_r', col_cluster=False)
heatmapPlotPdf = "{0}_heatmap_rearranged_ohne_TNF_Tgfb1.pdf".format(get_file_info(output_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

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

# Plot two heatmaps on the basis of overall means of transcription factors
# In [24]: peaksDF.drop('groups', axis=0).mean(axis=1).sort_values()                                                                                            
# name
# Il7        3.240490
# Il17a      3.560606
# Il5        3.975758
# IL6        4.161616
# Il25       4.636364
# Il2        4.774892
# Il10       4.949495
# Il21       5.318182
# Il22       5.479798
# Il23a      6.121212
# Il4        7.193182
# Ifng       8.186869
# Ztbt16     9.187328

# Il13      10.189394
# Tgfb1     13.063973
# RORc      14.350816
# TNF       21.151515

# Get the two dataframes for low and high values to plot in separate heatmaps
lowPeaksDF      = peaksDF.drop(['Tgfb1','RORc','TNF','groups']).astype('float')
highPeaksDF     = peaksDF.loc [['Tgfb1','RORc','TNF']].astype('float')
meanLowPeaksDF  = meanPeaksDF.drop(['Tgfb1','RORc','TNF']).astype('float')
meanHighPeaksDF = meanPeaksDF.loc [['Tgfb1','RORc','TNF']].astype('float')

# Plot heatmap for low values
g = sns.clustermap(lowPeaksDF, z_score=1, cmap='RdBu_r', col_cluster=False, vmin=-1.7, vmax=1.7, figsize=(20,8)); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
heatmapPlotPdf = "{0}_low_value_peaks_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')

g = sns.clustermap(highPeaksDF, z_score=1, cmap='RdBu_r', col_cluster=False, vmin=-1.7, vmax=1.7, figsize=(15,2)); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
heatmapPlotPdf = "{0}_high_value_peaks_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')

g = sns.clustermap(meanLowPeaksDF, z_score=1, cmap='RdBu_r', col_cluster=False, vmin=-1.7, vmax=1.7, figsize=(10,10)); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
heatmapPlotPdf = "{0}_mean_low_value_peaks_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')

g = sns.clustermap(meanHighPeaksDF, z_score=1, cmap='RdBu_r', col_cluster=False, vmin=-1.7, vmax=1.7, figsize=(10,3)); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
heatmapPlotPdf = "{0}_mean_high_value_peaks_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')


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

