# pwd
cd /home/rad/users/gaurav/users/thorsten/foxp1Loci

# Input file source: Original foxp1 chipseq publication (GSE31006)
# Parse fasta header to get foxp1 peaks location
# Take the peaks midpoint and get the start and end by -/+ 50 bp
# >P-146-1~chr10~123822928~ChIP:~170.682~control:~1.128~region:~123822692-123823121~ef:~151.313~ps:~65~cor:~0.593556~-log10_qv:~7632.53~-log10_pv:~7637.52~qv_rank:~33
# Extract chr10~123822928
# Chrom: chr10
# Start: 123822928 - 50
# End  : 123822928 + 50
# Name : P-146-1

# 1) ################################################################################################################
ipython

# Input file
input_file  = '/home/rad/users/gaurav/users/thorsten/foxp1Loci/input/GSM768310_FOXP1_ChIP-Seq_1_peaks_200bp.fa'

# Output file and handle
output_file = '/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/GSM768310_FOXP1_ChIP-Seq_1_peaks_200bp.bed'
ofile = open(output_file, 'w')

# Parse fasta file and generate the bed file
with open(input_file) as fp:
    for header, _seq in read_fasta(fp):
      name, chrom, midpoint, _rest = header.split('~', 3) # '>P-146-1', 'chr10', '123822928'
      start = int(midpoint) - 50 
      end   = int(midpoint) + 50 
      name  = "{0}_{1}_{2}_{3}".format(chrom, start, end, name.lstrip('>'))
      line  = "\t".join([chrom, str(start), str(end), name])
      ofile.write("{0}\n".format(line))
ofile.close()

# Sort the output bed file
os.system("sort -k1,1 -k2n {0} -o {0}".format(output_file))
print ("\n- Your text  output file is: {0}".format(output_file))

# Get the peaks tab file for atacseq data
# Input and output files
input_file = "/media/rad/SSD1/atac_temp/christine/AGRad_ATACseq_MUC001/analysis/normalized/5320_53631_53646_6075_merge_master_peaks_vstCounts.matrix"
peaksDF    = pd.read_csv(input_file, sep="\t", index_col=0)
# Save the topPeaks files in tab separated files
tabPeaksDF = peaksDF.copy()

# Get the chr, start and end columns
tabPeaksDF['tabDFSlice'] = tabPeaksDF.index.values.tolist()
tabDFSlice               = tabPeaksDF['tabDFSlice'].str.split("_", n = 3, expand = True)
tabPeaksDF['chr']        = tabDFSlice[0]
tabPeaksDF['start']      = tabDFSlice[1]
tabPeaksDF['end']        = tabDFSlice[2]

# Rearrange columns to new order
new_columns_order = ['chr', 'start', 'end', '5320_livmet-1_005_atac_030_s11', '5320_livmet-3_005_atac_031_s12', '5320_lungmet-1_005_atac_029_s10', '5320_ppt-1_004_atac_023_20k_s18', '5320_ppt-1_005_atac_028_s9', '53631_livmet-1_005_atac_025_s6', '53631_lungmet-2_005_atac_026_s7', '53631_lungmet-3_005_atac_027_s8', '53631_ppt-1_005_atac_024_s5', '53646_livmet-1_005_atac_033_s2', '53646_livmet-2_005_atac_034_s3', '53646_livmet-3_005_atac_035_s4', '53646_ppt-1_004_atac_022_s17', '53646_ppt-1_005_atac_032_s1', '53646_ppt-1_frozen_004_atac_022_fr_s16', '6075_blood_005_atac_038_s15', '6075_lungmet-1_005_atac_037_s14', '6075_ppt-1_005_atac_036_s13']
tabPeaksDF        = tabPeaksDF[new_columns_order]

# Save the DF as tab separated file
taboutput_file = "/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/5320_53631_53646_6075_merge_master_peaks_vstCounts.tab"
tabPeaksDF.to_csv(taboutput_file, sep='\t', index = False, float_format='%.2g')

############ USER DEFINIED CODE ############
def read_fasta(fp):
  ''' 
    - Parse a fasta file using a generator 
    - Source: https://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
  '''
  name, seq = None, []
  for line in fp:
      line = line.rstrip()
      if line.startswith(">"):
          if name: yield (name, ''.join(seq))
          name, seq = line, []
      else:
          seq.append(line)
  if name: yield (name, ''.join(seq))

# 2) ################################################################################################################
# Intersect with atacseq peaks

# 1) Get the peaks in the foxp1 region
peakMatrix="/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/5320_53631_53646_6075_merge_master_peaks_vstCounts.tab"

# 2.2) Get bed files for each sample
od="$(dirname ${peakMatrix})/samplesBed"
mkdir -p ${od}
header=$(head -1 ${peakMatrix})
for i in {4..21};
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

# 2.3) Get the mean score of atacseq peaks for the foxp1 region
pd="$(dirname ${peakMatrix})/foxp1Peaks"; mkdir -p ${pd}
ib="/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/GSM768310_FOXP1_ChIP-Seq_1_peaks_200bp.bed"
sort -k1,1 -k2n ${ib} -o ${ib}
for bbed in ${od}/*.bed;
do
 echo ${bbed}
 bedtools map -a ${ib} -b ${bbed} -c 4 -o mean > ${pd}/$(basename ${bbed} .bed)_foxp1_peaks_mean.bed
 # Add header to the file
 header="chr\tstart\tend\tname\t$(basename ${bbed} .bed)"
 sed -i "1i${header}" ${pd}/$(basename ${bbed} .bed)_foxp1_peaks_mean.bed
done

# Merge the individual mean peaks into a dataframe
ipython
from scipy.stats import zscore

# Read and merge into df
pkDF = pd.concat([pd.read_csv(f, sep='\t').set_index(['chr', 'start', 'end', 'name']) for f in glob.glob('/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/foxp1Peaks/*.bed')],axis=1).reset_index()
pkDF.replace('.',0, inplace=True)

# Plot the heatmap
peaksDF = pkDF.copy()
peaksDF.drop(columns = ['chr','start','end'], inplace = True)
peaksDF.set_index('name',inplace=True)
peaksDF = peaksDF.loc[~(peaksDF==0).all(axis=1)] # Remove rows with all 0s
peaksDF = peaksDF.astype('float') # Conver all columns to float
# Remove _peaks.bed_mean from column names
peaksDF = peaksDF.rename(columns={col: col.split('_peaks')[0] for col in peaksDF.columns})

# Drop the 004 samples
peaksDF.drop(list(peaksDF.filter(regex = '004')), axis = 1, inplace = True)

# Save to
output_file = '/home/rad/users/gaurav/users/thorsten/foxp1Loci/output/foxp1loci_merge_peaks_averageReads.txt'
peaksDF.to_csv(output_file, sep='\t', index = True, float_format='%.2g')

# Calculate the zscore of the dataframe
# The peaksDF.apply(zscore, axis=1) returns an array 
# https://stackoverflow.com/questions/35491274/pandas-split-column-of-lists-into-multiple-columns
from scipy.stats import zscore
peaksDF = peaksDF.astype('float')
zscorePeaksDF     = pd.DataFrame(pd.DataFrame(peaksDF.apply(zscore, axis=1))[0].values.tolist(), columns=peaksDF.columns.tolist(), index=peaksDF.index)

# Plot the heatmap for all replicates separately for all time points
sns.set(font_scale=0.25)
g = sns.clustermap(zscorePeaksDF, cmap='RdBu_r', figsize=(5,15)); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0); cax = plt.gcf().axes[-1]; cax.tick_params(labelsize=5);
heatmapPlotPdf = "{0}_heatmap.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')

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
