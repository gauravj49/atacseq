#!/usr/local/bin/python	
"""	
***********************************************	
- PROGRAM: customFeatureHeatmap.py	
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)	
***********************************************	
"""	

def main():	
  print(os.system("python -VV"))
  # feature_name     = "ETP_Class_small"
  # input_dir = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/cisPeaks"
  # output_file = "/media/rad/HDD1/atacseq/anja/tALLcellLineMm/analysis/mergedReps/cis/{}Reads.txt".format(feature_name)

  # Read and merge into df
  pkDF = pd.concat([pd.read_csv(f, sep='\t').set_index(['chr', 'start', 'end','name']) for f in glob.glob("{0}/*{1}_peaks_mean.bed".format(input_dir, feature_name))],axis=1).reset_index()
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
  parser.add_argument("-id", metavar="--indir" , help="*Input directory containing the counts(peaks) files", dest="input_dir"  , type=str, required=True)
  parser.add_argument("-sf", metavar="--cfile" , help="*input sample files list of counts(peaks) files separated by comma. Example: /path/a.txt,/path/b.txt", dest="sample_names_file" , type=str, required=True)
  parser.add_argument("-of", metavar="--ofile" , help="*Output file to save the annotation and plots file", dest="output_file", type=str, required=True)
  parser.add_argument("-ft", metavar='--ftname', help="*Name of the feature. Ex. ETPsmall or CTCFsites", dest="feature_name", type=str, required=True)

  parser.add_argument("-ff", metavar="--ffile" , help=" List of selected element one per line. Ex: gene names or cis sites" , dest="element_file", type=str)
  parser.add_argument("-fl", metavar="--ftlist", help=" Comma separated list of selected elements to plot\n - (Default: [])", dest="selectedElementList" , type=str, default = list())
  parser.add_argument("-ds", metavar="--dropsm", help=" Comma separated list of samples to drop.\n - (Default: [])" , dest="dropSamples" , type=str, default = list())
  parser.add_argument("-mn", metavar="--vmin"  , help=" Set the min value of the heatmap colorbar. Ex: -2.5" , dest="vmin", type=str, default="")
  parser.add_argument("-mx", metavar="--vmax"  , help=" Set the max value of the heatmap colorbar. Ex:  2.5" , dest="vmax", type=str, default="")
  parser.add_argument('-xc', "--xclustering"   , help="if set, cluster x-rows", action='store_true', default=False)
  parser.add_argument('-yc', "--yclustering"   , help="if set, cluster y-rows", action='store_true', default=False)

  # Print the help message only if no arguments are supplied	
  if len(sys.argv)==1:	
      parser.print_help()	
      sys.exit(1)	

  # Save the STDOUT output in a log file	
  if parser.parse_args().input_dir:	
      logdir="{0}/logs".format(parser.parse_args().input_dir)	
      create_dir(logdir)	
      logfile = "{0}/{1}.log".format(logdir, parser.parse_args().feature_name)	
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

  # Get input parameters
  input_dir         = args.input_dir.rstrip('\/')
  sample_names_file = args.sample_names_file
  output_file       = args.output_file
  xclustering       = args.xclustering
  yclustering       = args.yclustering
  vmin              = args.vmin
  vmax              = args.vmax
  selected_features = args.selected_features   # Selected features like clinical features to keep
  selected_featfile = args.selected_featfile   # Selected features in a file ... one feature per line 
  if selected_features:
      selected_features = selected_features.split(",")
      if verbose:
          print "Selected_features: "
          for s in selected_features:
              print "\t- {0}".format(s)

  selected_features = list()
  if selected_featfile:
        with open(selected_featfile,'rU') as fg:
            # head
            # ENSMUSG00000025130
      # ENSMUSG00000025393
            # ....
            selected_features = [line for line in useful_lines(fg)]

  # Create output directory if not present
  create_dir(get_file_info(output_file)[0])

  # User variables
  # input_dir = "{0}/cisPeaks".format(input_dir)
  # output_file = "{0}/{1}_peaks_averageReads.txt".format(input_dir, feature_name)

  main()
