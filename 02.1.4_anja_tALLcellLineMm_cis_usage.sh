# PWD
cd /home/rad/users/gaurav/projects/seqAnalysis/atacseq

# 1) Get the peaks in the cis region
peakMatrix="/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/tALLcellLineMm_mergedreps_all_merge_master_peaks_vstCounts.tab"
cisDir=$(dirname ${peakMatrix})/cis

# 2.2) Get bed files for each sample
od="${cisDir}/samplesBed"
mkdir -p ${od}
header=$(head -1 ${peakMatrix})
numCols=$(awk -F'\t' '{print NF; exit}' ${peakMatrix})
for i in $(eval echo "{4..${numCols}}")
do
 bname=$(cut -f${i} <<< ${header})
 peaksBedFile="${od}/${bname}_peaks.bed"
 echo ${bname}
 cut -f1-3,${i} ${peakMatrix} > ${peaksBedFile}

 # Remove the header
 sed -i '1d' ${peaksBedFile}

 # Sort the input file with the header line intact
 # (head -n 1 ${peaksBedFile} && tail -n +2 ${peaksBedFile} | sort -k1,1V -k2,2g -k3,3g) > ${peaksBedFile}.tmp && mv ${peaksBedFile}.tmp ${peaksBedFile}
 sort -k1,1V -k2,2g -k3,3g ${peaksBedFile} -o ${peaksBedFile}
done

# Remove random chromosome GL
rm -rf ${od}/PeakID_peaks.bed
for f in ${od}/*.bed; do sed -i '/GL/d' ${f}; done;
for f in ${od}/*.bed; do sed -i '/JH/d' ${f}; done;

# 2.3) Get the mean score of atacseq peaks for the cis region
pd="${cisDir}/cisPeaks"; mkdir -p ${pd}
ib="/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/ETP_Class_small.bed"
ie="/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/ETPsmall.bed"
ic="/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/classSmall.bed"
gf="/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/mmuGenome.info"
sort -k1,1V -k2,2g -k3,3g ${ib} -o ${ib}

# Remove the outlier cis
grep -v "ETPsmall_chr11_3183653_3194348" ${ib} > tmp && mv tmp ${ib}
grep -v "ETPsmall_chr11_3183653_3194348" ${ie} > tmp && mv tmp ${ie}
for bbed in ${od}/*.bed;
do
 echo ${bbed}
 bedtools map -a ${ib} -b ${bbed} -c 4 -o mean -g ${gf} > ${pd}/$(basename ${bbed} .bed)_ETP_Class_small_peaks_mean.bed
 bedtools map -a ${ic} -b ${bbed} -c 4 -o mean -g ${gf} > ${pd}/$(basename ${bbed} .bed)_classSmall_peaks_mean.bed
 bedtools map -a ${ie} -b ${bbed} -c 4 -o mean -g ${gf} > ${pd}/$(basename ${bbed} .bed)_ETPsmall_peaks_mean.bed
 # Add header to the file
 header="chr\tstart\tend\tname\t$(basename ${bbed} _mmu_atacseq_pe_R1_rmdup_peaks.bed)_mean"
 sed -i "1i${header}" ${pd}/$(basename ${bbed} .bed)_ETP_Class_small_peaks_mean.bed
 sed -i "1i${header}" ${pd}/$(basename ${bbed} .bed)_classSmall_peaks_mean.bed
 sed -i "1i${header}" ${pd}/$(basename ${bbed} .bed)_ETPsmall_peaks_mean.bed
done

# Plot heatmaps
python scripts/customFeatureHeatmap.py -cd=/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis -ct=ETPsmall
python scripts/customFeatureHeatmap.py -cd=/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis -ct=classSmall
python scripts/customFeatureHeatmap.py -cd=/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis -ct=ETP_Class_small

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

