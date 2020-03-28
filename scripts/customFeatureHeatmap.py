#!/usr/local/bin/python	
"""	
***********************************************	
- PROGRAM: annotate_variant_bash_wrapper.py	
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)	
***********************************************	
"""	

def main():	
  # cisType     = "ETP_Class_small"
  # cisPeaksDir = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/cisPeaks"
  # output_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/{}Reads.txt".format(cisType)

  # # Create a meta log file
  # metaLogFile = "{0}/logs/{1}_meta.log".format(cisDir, cisType); create_dir("{0}/logs".format(cisDir))

  # # Log everything in the meta log file
  # # https://serverfault.com/questions/103501/how-can-i-fully-log-all-bash-scripts-actions
  # # Explanation:
  # # 1.  exec 3>&1 4>&2
  # #     Saves file descriptors so they can be restored to whatever they were before redirection 
  # #     or used themselves to output to whatever they were before the following redirect.
  # # 2.  trap 'exec 2>&4 1>&3' 0 1 2 3
  # #     Restore file descriptors for particular signals. Not generally necessary since they should 
  # #     be restored when the sub-shell exits.
  # # 3.  exec 1>log.out 2>&1
  # #     Redirect stdout to file log.out then redirect stderr to stdout. Note that the order is 
  # #     important when you want them going to the same file. stdout must be redirected before stderr is redirected to stdout.

  # ofile.write("#!/bin/bash\n\nexec 3>&1 4>&2\ntrap 'exec 2>&4 1>&3' 0 1 2 3\nexec 1> {} 2>&1\n".format(metaLogFile))
  # # Everything below will go to the log file

  # pyVersion = os.system("python -VV")
  # ofile.write("\n########################################\n# Log file: logs/{0}_meta.log\n{1}\n########################################\n".format(sampleName, pyVersion))
  print(os.system("python -VV"))

  # cisType     = "ETP_Class_small"
  # cisPeaksDir = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/cisPeaks"
  # output_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/{}Reads.txt".format(cisType)

  # Read and merge into df
  pkDF = pd.concat([pd.read_csv(f, sep='\t').set_index(['chr', 'start', 'end','name']) for f in glob.glob("{0}/*{1}_peaks_mean.bed".format(cisPeaksDir, cisType))],axis=1).reset_index()
  pkDF.replace('.',0, inplace=True)

  # Save to output file
  pkDF.to_csv(output_file, sep='\t', index = False, float_format='%.2g')

  # Plot the heatmap
  peaksDF = pkDF.copy()
  peaksDF.drop(columns = ['chr','start','end'], inplace = True)
  peaksDF.set_index('name',inplace=True)
  peaksDF = peaksDF.loc[~(peaksDF==0).all(axis=1)] # Remove rows with all 0s
  peaksDF = peaksDF.astype('float') # Conver all columns to float
  peaksDF = peaksDF.rename(columns={col: col.split('_peaks')[0] for col in peaksDF.columns})

  # Rearrange peaksDF columns
  new_columns_order = ['TransPB_20191011114994_P00000011_HSC-BoneMarrow_mean','TransPB_20191011114931_P00000011_MPP-BoneMarrow_mean','TransPB_20191011108131_P00000011_CLP-BoneMarrow_mean','TransPB_20191011114755_P00000011_ETP-Tcell_mean','TransPB_20191011114671_P00000011_DN2a-Tcell_mean','TransPB_20191011143464_P00000011_DN2b-Tcell_mean','TransPB_20191011072548_P00000011_DN3-Tcell_mean','TransPB_20191011090063_P00000011_DN4-Tcell_mean','TransPB_20191011125071_P00000011_DP-Tcell_mean','TransPB_20191011142220_P00000011_CD4-SinglePositive_mean','TransPB_20191011143412_P00000011_CD8-SinglePositive_mean','TransPB_20191011137800_P00000011_CD8-BoneMarrow_mean','TransPB_20191011141045_P00000011_Bcell-BoneMarrow_mean','TransPB_20191011097119_P00000011_NK-BoneMarrow_mean']

  groups = ['HSC','MPP','CLP','ETP','DN2a','DN2b','DN3','DN4','DP','CD4','CD8SP','CD8BM','Bcell','NK']

  # Rearrange columns to new order
  peaksDF = peaksDF[new_columns_order]
  peaksDF = peaksDF.astype('float') # Conver all columns to floa

  # Make column names smaller
  peaksDF.columns = groups

  # Plot the dataframe as heatmap
  # yticklabels=1 to show every ticklabel
  fig = plt.figure(figsize=(20,15))
  g   = sns.clustermap(peaksDF, z_score=0, cmap='RdBu_r', col_cluster=False, square=True, yticklabels=1)
  g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 3)
  g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 3)
  heatmapPlotPdf = "{0}_heatmap.pdf".format(get_file_info(output_file)[3])
  plt.tight_layout()
  plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
  plt.close('all')
  print(heatmapPlotPdf)

  # # Calculate the zscore of the dataframe
  # # The peaksDF.apply(zscore, axis=1) returns an array 
  # # https://stackoverflow.com/questions/35491274/pandas-split-column-of-lists-into-multiple-columns
  # from scipy.stats import zscore
  # peaksDF = peaksDF.astype('float')
  # zscorePeaksDF     = pd.DataFrame(pd.DataFrame(peaksDF.apply(zscore, axis=1))[0].values.tolist(), columns=peaksDF.columns.tolist(), index=peaksDF.index)

  # # Plot the heatmap for all replicates separately for all time points
  # g = sns.clustermap(zscorePeaksDF, cmap='RdBu_r', col_cluster=False, figsize=(20,15), square=True, yticklabels=1); plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0); cax = plt.gcf().axes[-1]; cax.tick_params(labelsize=5);
  # heatmapPlotPdf = "{0}_heatmap_zscore.pdf".format(get_file_info(output_file)[3]); plt.savefig(heatmapPlotPdf,bbox_inches = 'tight'); plt.close('all')
  # print(heatmapPlotPdf)

################ USER DEFINED FUNCTIONS ###################	

def check_options():	
  ''' Checks the options to the program '''	

  # Create parser object	
  parser = argparse.ArgumentParser(add_help=True,	
      formatter_class=argparse.RawTextHelpFormatter,	
      epilog=textwrap.dedent('''
      ----------------- SAMPLE USAGE ------------------	
      python scripts/customFeatureHeatmap.py -cd=/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis -ct=ETPsmall	
      -------------------------------------------------	
      CONTACT: 	
          Gaurav Jain	
          gaurav.jain@tum.de	
      -------------------------------------------------	
      '''))	

  # Add arguments 	
  parser.add_argument("-cd", metavar='--cisdir', help="cis directory", dest="cisDir" , type=str, required=True)	
  parser.add_argument("-ct", metavar='--cistyp', help="type of cis"  , dest="cisType", type=str, required=True)	
  

  # http://thomas-cokelaer.info/blog/2014/03/python-argparse-issues-with-the-help-argument-typeerror-o-format-a-number-is-required-not-dict/	

  # Print the help message only if no arguments are supplied	
  if len(sys.argv)==1:	
      parser.print_help()	
      sys.exit(1)	

  # Save the STDOUT output in a log file	
  if parser.parse_args().cisDir:	
      logdir="{0}/logs".format(parser.parse_args().cisDir)	
      create_dir(logdir)	
      logfile = "{0}/{1}.log".format(logdir, parser.parse_args().cisType)	
  else:	
      logdir  = "{0}/logs".format(os.getcwd())	
      create_dir(logdir)	
      logfile = "{0}/{1}.log".format(logdir,get_file_info(sys.argv[0])[1])	

  logf = open(logfile, 'w')	
  sys.stdout = Log(logf, sys.stdout)	

  # Parse command line with parse_args and store it in an object	
  args = parser.parse_args()	
  print_initial_arguments(parser)	
  return args	

if __name__=="__main__":	
  from functools import partial, reduce	
  print (__doc__)	

  # Built in modules	
  import argparse	
  import os.path	
  import sys	

  # 3rd party modules	
  import textwrap	
  import re	
  import numpy as np	
  import scipy as sp	
  import pandas as pd	
  import matplotlib as mp	
  #mp.use('Agg') # to use matplotlib without X11	
  import matplotlib.pyplot as plt	
  import subprocess	
  import binascii as bi	
  import scipy.stats as stats	
  from collections import *	
  from numpy import nanmean
  from scipy.stats import zscore

  # for looping files in a dir	
  import glob	

  # user defined modules	
  from gjainPyLib import *     # import all the functions from the Gaurav`s python scripts/library	

  ### for color scale	
  from  matplotlib import colors	
  from itertools import cycle, islice # barplot colors	

  ################ USER CONFIGURATION ###################	
  np.set_printoptions(precision=6)	
  #######################################################	

  # Get input options	
  args = check_options()	

  # Store the variables
  cisType = args.cisType
  cisDir  = args.cisDir

  # User variables
  cisPeaksDir = "{0}/cisPeaks".format(cisDir)
  output_file = "{0}/{1}_peaks_averageReads.txt".format(cisDir, cisType)

  main()
