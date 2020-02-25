#!/usr/local/bin/python	
"""	
***********************************************	
- PROGRAM: 00_run_atacseq_wrapper.py	
- CONTACT: Gaurav Jain(gaurav.jain@tum.de)	
***********************************************	
"""	

def main1():	
    pass

def main1():	

    # Project folder structure
    # ├── analysis
    # ├── annotation
    # ├── data
    # │   ├── bams
    # │   ├── bed
    # │   ├── bigwigs
    # │   └── fastq
    # ├── logs
    # ├── peaks
    # └── qc

    # Get output file name	
    scriptDir  = "{0}/scripts/wrapper/{1}".format(os.getcwd(), projName); os.system("mkdir -p {0}".format(scriptDir))	
    scriptFile = "{0}/{1}_atacseq_wrapper.sh".format(scriptDir, sampleName)
    ofile      = open(scriptFile, 'w')	

    # Create a meta log file
    metaLogFile = "{0}/logs/{1}_atacseq_pipeline.log".format(projDir, sampleName); create_dir("{0}/logs".format(projDir))

    # Log everything in the meta log file
    # https://serverfault.com/questions/103501/how-can-i-fully-log-all-bash-scripts-actions
    # Explanation:
    # 1.  exec 3>&1 4>&2
    #     Saves file descriptors so they can be restored to whatever they were before redirection 
    #     or used themselves to output to whatever they were before the following redirect.
    # 2.  trap 'exec 2>&4 1>&3' 0 1 2 3
    #     Restore file descriptors for particular signals. Not generally necessary since they should 
    #     be restored when the sub-shell exits.
    # 3.  exec 1>log.out 2>&1
    #     Redirect stdout to file log.out then redirect stderr to stdout. Note that the order is 
    #     important when you want them going to the same file. stdout must be redirected before stderr is redirected to stdout.

    ofile.write("#!/bin/bash\n\nexec 3>&1 4>&2\ntrap 'exec 2>&4 1>&3' 0 1 2 3\nexec 1> {} 2>&1\n".format(metaLogFile))
    # Everything below will go to the log file

    ofile.write("\n########################################\n# Log file: logs/{0}_atacseq_pipeline.log\npython -VV\n########################################\n".format(sampleName))

    # 1) Pre-analysis (quality control (QC) and alignment)
    ofile.write("\n#- 1.1) Convert merged vcf file to bed file\n")
    # bash scripts/01_map_singleendFastQ_bowtie2.sh ${fastqdir} ${outputDir} ${projName} ${species}









    # 1.1) Convert merged vcf file to bed file	
    ofile.write("\n#- 1.1) Convert merged vcf file to bed file\n")	
    snvBedFile     = "{}.bed".format(interimFileName)	
    snvUcscBedFile = "{}.ucscbed".format(interimFileName)	
    ofile.write("time python scripts/export_variant_input_to_bed.py -if={0} -of={1}\n".format(inputVariantFile, interimFileName))	
    ofile.write("egrep -v 'Un|random|GL|JH|M|KQ|KZ' {0} > {1} && mv {1} {0}\n".format(snvBedFile,"{0}.tmp".format(snvBedFile))) 	
    ofile.write("egrep -v 'Un|random|GL|JH|M|KQ|KZ' {0} > {1} && mv {1} {0}\n".format(snvUcscBedFile,"{0}.tmp".format(snvUcscBedFile)))	
    ofile.write("rm -rf {0} {1}\n".format("{0}.tmp".format(snvBedFile), "{0}.tmp".format(snvUcscBedFile)))	

    # 1.2) Get the the subset bam for variant reads only
    ofile.write("\n#\t- 1.2) Get the the subset bam for variant reads only\n")	
    snvSubBamFile = "{}_subset_variant.bam".format(interimFileName)	
    ofile.write("time intersectBed -abam {0} -b {1} -wa > {2}\n".format(tumorBamFile, snvBedFile, snvSubBamFile))	

    # 1.3) Make genome information (chromosome sizes) from the BAM file to get a genome.info file	
    genomeInfoFile     = "{}_genome.info".format(interimFileName)	
    genomeInfoUcscFile = "{}_genome.ucscinfo".format(interimFileName)	
    ofile.write("\n#- 1.3) Make genome information (chromosome sizes) from the BAM file to get a genome.info file\n")	
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print $1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ' | sort -k1,1V -k2,2g -k3,3g > {1}\n".format(tumorBamFile, genomeInfoFile))	
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print \"chr\",$1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ' | sort -k1,1V -k2,2g -k3,3g > {1}\n".format(tumorBamFile, genomeInfoUcscFile))	

    # 2) Get the base reads quality and mapping quality for the variants	
    ofile.write("\n#- 2) Get the base reads quality and mapping quality for the variants\n")	
    # Get the pileup information for the variants 	
    ofile.write("\n#\t- Get the pileup information for the variants\n")	
    snvMpileupFile  =  "{}.pileup".format(interimFileName)	
    ofile.write("time samtools mpileup -A -B -Q 0 -s -x --output-QNAME -l {0} -f {1} {2} -o {3}".format(snvBedFile, genomeFastaFile, tumorBamFile,snvMpileupFile))	
    # 2.2) Parse the output of pileup	
    snvPileupFile  = "{}_parsed_tumor_exact_pileup.txt".format(interimFileName)	
    ofile.write("\n#\t- Get the pileup information for the variants\n")	
    ofile.write("time python scripts/parse_mpileup.py -ip={0} -ib={1}\n".format(snvMpileupFile, snvBedFile))	

    # 3) Segmental duplication	
    ofile.write("\n#- 3) Get the segmental duplication information\n")	
    segmentalDupBedFile = "{}_segmentalDups.bed".format(interimFileName)	
    tempSegDupfile      = get_temp_file()	
    ofile.write("egrep -v 'Un|random|GL|JH|M|KQ|KZ' {0} > {1}\n".format(genomeSegDupFile, tempSegDupfile))	
    ofile.write("time intersectBed -a {0} -b {1} -c -g {3}| sed 's/chr//' > {2}\n".format(snvUcscBedFile, tempSegDupfile, segmentalDupBedFile, genomeInfoUcscFile))	
    ofile.write("\n#\t- Adding the header to the file\n")	
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tSegmentalDupsCount' {}\n".format(segmentalDupBedFile))	

    # 4) Repeat masker annotation	
    ofile.write("\n#- 4) Get the repeat masker annotation\n")	
    repeatMaskerBedFile = "{}_repeatmasker.txt".format(interimFileName)	
    ofile.write("time intersectBed -a {0} -b {1} -wao -g {3} > {2}\n".format(snvBedFile, genomeRepMaskFile, repeatMaskerBedFile, genomeInfoFile))	
    ofile.write("\n#\t- Adding the header to the file\n")	
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tRepChrom\\tRepGenoStart\\tRepGenoEnd\\tRepNameClassFamily\\tSwScore\\tSwStrand\\tRepMaskFlag' {}\n".format(repeatMaskerBedFile))	

    # 5) GC percent annotation	
    ofile.write("\n#- 5) GC percent annotation\n")	
    # 5.1) For snvBinGCpc	
    ofile.write("\n#\t- Get the GC% for the 5bp bin the SNV overlaps\n")	
    snvBinGCpcFile = "{}_5bp_GC_percentage.txt".format(interimFileName)	
    ofile.write("time intersectBed -a {0} -b {1} -wb -g {3} > {2}\n".format(genome5bpGCpcFile, snvUcscBedFile, snvBinGCpcFile, genomeInfoUcscFile))	
    ofile.write("sed -i 's/chr//g' {}".format(snvBinGCpcFile))	
    ofile.write("\n#\t- Adding the header to the file\n")	
    ofile.write("sed -i '1iGCpc5bpChrom\\tGCpc5bpStart\\tGCpc5bpEnd\\tGCVariantBasePercentage\\tChrom\\tStart\\tGenomicPos\\tVariantIDKey' {}\n".format(snvBinGCpcFile))	

    # 5.2) For snvReadMeanGCpc	
    ofile.write("\n#\t- Get the mean GC% for all the reads the SNV overlaps\n")	
    ofile.write("\n#\t- Map the reads to the snv locations, concatenate all the reads and calculate the mean gc% of all the mapped reads\n")	
    snvReadMeanGCpcFile = "{}_snvReads_Mean_GC_percentage.txt".format(interimFileName)	

    ofile.write('time bedtools map -prec 2 -a {0} -b {1} -c 10,10 -o count,concat -g {3} | awk -v OFS="\\t" \'{{n=length($6); gc=gsub("[gcGC]", "", $6); print $1,$2,$3,$4,gc/n}}\'  > {2}\n'.format(snvBedFile, snvSubBamFile, snvReadMeanGCpcFile, genomeInfoFile))	
    ofile.write("\n#\t- Adding the header to the file\n")	
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tVariantReadsMeanGCPc' {}\n".format(snvReadMeanGCpcFile))	

    # 6) Get the 10 bp sequence around the SNV from the genome	
    ofile.write("\n#- 6) Get the 10 bp sequence around the SNV from the genome\n")	
    tenBpAroundSNVFile  = "{}_10bp_around_SNV.txt".format(interimFileName)	
    temp10bpbinBaseFile = get_temp_file() 	
    ofile.write("time awk -v OFS=\"\\t\" '{{print $1,$2-10,$3+10,$3}}' {0} > {1}\n".format(snvBedFile, temp10bpbinBaseFile))	
    ofile.write("time bedtools getfasta -fi {0}  -bed {1} -name -bedOut > {2}\n".format(genomeFastaFile, temp10bpbinBaseFile, tenBpAroundSNVFile))	
    ofile.write("\n#\t- Adding the header to the file\n")	
    ofile.write("sed -i '1iChrom\\tStart\\tEnd\\tGenomicPos\\tVariant10bpFlankingSequence' {}\n".format(tenBpAroundSNVFile))	

    # 7) Get the distance to the closest probe from the SNV
    # WESProbeAnnotation = wesProbe_probeID|TargetID|Sequence|Replication
    ofile.write("\n#- 7) Get the distance to the closest probe from the SNV\n")	
    closestProbeFile = "{}_closest_wes_probe.txt".format(interimFileName)	
    tempclosestprobe = get_temp_file()
    tmpsortingfile   = get_temp_file()
    ofile.write("grep -v ^M {0} > {1}\n".format(wesProbesFile, tempclosestprobe))	
    ofile.write("(head -n 2 {0} && tail -n +3 {0} | sort -k1,1V -k2,2g -k3,3g) > {1} && mv {1} {0}\n".format(tempclosestprobe, tmpsortingfile))	
    ofile.write("time bedtools closest -a {0} -b {1} -d -t first -g {3}| uniq > {2}\n".format(snvBedFile, tempclosestprobe, closestProbeFile, genomeInfoFile))	
    ofile.write("\n#\t- Adding the header to the file\n")	
    ofile.write("sed -i '1iChrom\\tStart\\tGenomicPos\\tVariantIDKey\\tWESProbeChrom\\tWESProbeStart\\tWESProbeEnd\\tWESProbeAnnotation\\tWESProbeScore\\tWESProbeStrand\\tWESProbeDistanceToVariant' {}\n".format(closestProbeFile))	

    # 8.1) Make genome information (chromosome sizes) from the Normal BAM file to get a genome.info file	
    normGenomeInfoFile     = "{}_genome_normal.info".format(interimFileName)	
    normGenomeInfoUcscFile = "{}_genome_normal.ucscinfo".format(interimFileName)	
    ofile.write("\n#- 8.1) Make genome information (chromosome sizes) from the Normal BAM file to get a genome.info file\n")	
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print $1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ' > {1}\n".format(normalBamFile, normGenomeInfoFile))	
    ofile.write("samtools view -H {0}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print \"chr\",$1,\"\\t\",$2,\"\\n\"}}' | egrep -v 'Un|random|GL|JH|M|KQ|KZ' > {1}\n".format(normalBamFile, normGenomeInfoUcscFile))	

    # 8.2) Get usable indels list	
    ofile.write("\n#- 8.2) Get usable tumor indels information\n")	
    tumorIndelFile     = "{}_indels.txt".format(interimFileName)	
    usableIndelBed     = "{}_indels.bed".format(interimFileName)	
    usableIndelUcscBed = "{}_indels.ucscbed".format(interimFileName)	
    ofile.write("\n#- python scripts/get_usable_indels_for_annoatation.py -if={0} -pf={1} -bm={2} -od={3}\n".format(inputVariantFile, snvPileupFile, snvSubBamFile, interimSampleDir))	
    ofile.write("time python scripts/get_usable_indels_for_annoatation.py -if={0} -pf={1} -bm={2} -od={3}\n".format(inputVariantFile, snvPileupFile, snvSubBamFile, interimSampleDir))	

    # 8.3) Increase the size of each feature in 5 bases on both sizes	
    usableIndelExtendedBed     = "{}_indels.extendedbed".format(interimFileName)	
    usableIndelExtendedUcscBed = "{}_indels.extendeducscbed".format(interimFileName)	
    ofile.write("\n#- 8.3) Increase the size of each feature in 5 bases on both sizes\n")	
    ofile.write("time bedtools slop -i {0} -g {1} -b 5 > {2}\n".format(usableIndelBed, normGenomeInfoFile, usableIndelExtendedBed))	
    ofile.write("time bedtools slop -i {0} -g {1} -b 5 > {2}\n".format(usableIndelUcscBed, normGenomeInfoUcscFile, usableIndelExtendedUcscBed))	

    # 9.1) Get the base reads and mapping quality for the variants for indels in tumor bam files	
    ofile.write("\n#- 9.1) Get the base reads and mapping quality for the variants for indels in tumor bam files\n")	
    # Get the pileup information for the variants
    ofile.write("\n#\t- Get the pileup information for the variants in tumor bam file\n")	
    snvMextendedTumorPileupFile  =  "{}_tumor.extendedpileup".format(interimFileName)	
    ofile.write("time samtools mpileup -A -B -Q 0 -s -x --output-QNAME -l {0} -f {1} {2} -o {3}".format(usableIndelExtendedBed, genomeFastaFile, tumorBamFile,snvMextendedTumorPileupFile))

    # 9.2) Parse the output of extended pileup file for tumor bams	
    ofile.write("\n#- 9.2) Parse the output of extended pileup file for tumor bams\n")	
    extParsedTumorPileupFile  = "{}_tumor_parsed_extendedpileup.txt".format(interimFileName)	
    ofile.write("\n#\t- Get the extended pileup information for the variants\n")	
    ofile.write("time python scripts/parse_extended_mpileup.py -ip={0} -ib={1} -pt=t \n".format(snvMextendedTumorPileupFile, usableIndelBed))	

    # 9.3) Get the base reads and mapping quality for the variants for indels in normal bam files	
    ofile.write("\n#- 9.3) Get the base reads and mapping quality for the variants for indels in normal bam files")	
    # Get the pileup information for the variants 	
    ofile.write("\n#\t- Get the pileup information for the variants in normal bam file\n")	
    snvMextendedNormalPileupFile  =  "{}_normal.extendedpileup".format(interimFileName)	
    ofile.write("time samtools mpileup -A -B -Q 0 -s -x --output-QNAME -l {0} -f {1} {2} -o {3}".format(usableIndelExtendedBed, genomeFastaFile, normalBamFile,snvMextendedNormalPileupFile))	

    # 9.4) Parse the output of extended pileup file for normal bams	
    ofile.write("\n#- 9.4) Parse the output of extended pileup file for normal bams")	
    extParsedNormalPileupFile  = "{}_normal_parsed_extendedpileup.txt".format(interimFileName)	
    ofile.write("\n#\t- Get the extended pileup information for the variants\n")	
    ofile.write("time python scripts/parse_extended_mpileup.py -ip={0} -ib={1} -pt=n \n".format(snvMextendedNormalPileupFile, usableIndelBed))	

    # 10) Merge nnotations to the original file (mergedvcf)	
    ofile.write("\n#- 10) Merge annotations to the original file (mergedvcf)\n")	
    ofile.write("\n#\t- python scripts/merge_annotation_files.py -if={0} -pf={1} -sf={2} -rf={3} -gb={4} -ab={5} -gr={6} -cp={7} -ti={8} -tx={9} -nx={10}\n".format(inputVariantFile, snvPileupFile, segmentalDupBedFile, repeatMaskerBedFile, snvBinGCpcFile, tenBpAroundSNVFile, snvReadMeanGCpcFile, closestProbeFile, tumorIndelFile, extParsedTumorPileupFile, extParsedNormalPileupFile))	
    ofile.write("time python scripts/merge_annotation_files.py -if={0} -pf={1} -sf={2} -rf={3} -gb={4} -ab={5} -gr={6} -cp={7} -ti={8} -tx={9} -nx={10}\n".format(inputVariantFile, snvPileupFile, segmentalDupBedFile, repeatMaskerBedFile, snvBinGCpcFile, tenBpAroundSNVFile, snvReadMeanGCpcFile, closestProbeFile, tumorIndelFile, extParsedTumorPileupFile, extParsedNormalPileupFile))	
    ofile.close()	
    print("\n- Your text  output file is: {0}".format(scriptFile))	

################ USER DEFINED FUNCTIONS ###################	

def check_options():	
    ''' Checks the options to the program '''	

    # Create parser object	
    parser = argparse.ArgumentParser(
        prog='scripts/00_run_atacseq_wrapper.py', 
        usage='%(prog)s [options]',
        add_help=True,	
        formatter_class=argparse.RawTextHelpFormatter,	
        epilog=textwrap.dedent('''\	
        ----------------- SAMPLE USAGE ------------------	
        python scripts/00_run_atacseq_wrapper.py -h 	
        -------------------------------------------------	
        CONTACT: 	
                Gaurav Jain	
                gaurav.jain@tum.de	
        -------------------------------------------------	
        '''))	

    # Add arguments 	

    # Add the arguments for required ChIPseq arguments groups
    reqAtac_args = parser.add_argument_group('Required ATACseq Arguments')
    reqAtac_args.add_argument("-pd", metavar='--prjdir'    ,dest="projDir"   , type=str, help="* Project directory for output files and folders",  required=True)	
    reqAtac_args.add_argument("-ex", metavar="--expname"   ,dest="expName"   , type=str, help="* Experiment name to be used during the analysis.\n\t- Example: nkTimecourse, foxp1Loci etc.", required=True)
    reqAtac_args.add_argument("-sn", metavar='--sample'    ,dest="sampleName", type=str, help="* Name of the sample. ex: CD4PosTcells_mm_atacseq_se_R1" , required=True)
    reqAtac_args.add_argument("-sp", metavar="--species"   ,dest="species"   , type=str, help="* Species (organism) of the sample.\n\t- Supported species: {'hsa', 'mmu'}",  required=True, choices=['hsa','mmu'])


    reqAtac_args.add_argument("-tl", metavar="--trtList"   , help="* Text file containing ChIP sample filenames with full path (one samples per line)\n\t- Supported formats: fastq, SAM and BAM.", dest="treatmentList", type=str, required=True)
    reqAtac_args.add_argument("-pc", metavar="--peakCaller", help="* The peak calling software. \n\t- Choices: {'macs2narrow', 'macs2broad', 'bcp', 'epic'}", dest="peakCaller", type=str, required=True, choices=['macs2narrow','macs2broad','bcp','epic2'])

    # Add the arguments for optional ChIPseq groups
    optAtac_args = parser.add_argument_group('Optional ATACseq Arguments')
    optAtac_args.add_argument("-il", metavar="--inpList", help="Same as treatmentList but for the input/control files." , dest="inputList", default='', type=str)
    optAtac_args.add_argument("-qv", metavar="--qValue" , help="q-Value for peak calling.", dest="qValue", default='', type=str)
    optAtac_args.add_argument("-pv", metavar="--pValue" , help="p-Value for peak calling.", dest="pValue", default='', type=str)
    optAtac_args.add_argument("-du", metavar="--mrDup"  , help="Mark or remove duplicates.\n- Choices: {'mark', 'remove', 'none'}\n- Default: %(default)s", dest="mrDup" , default="mark" , type=str, choices=['mark','remove','none'])
    optAtac_args.add_argument("-bw", metavar="--bigWig" , help="BigWig creation for data visualization.\n- Choices: {'yes','no'}\n- Default: %(default)s", dest="bigWig", default="no", type=str, choices=['no','yes'])
    optAtac_args.add_argument("-bs", metavar="--binSize", help="Bin size for Bigwig creation (for resolution).\n- Default: %(default)s", dest="binSize", default="20", type=int)
    optAtac_args.add_argument("-nb", metavar="--normBW" , help="Normalization method for Bigwig creation.", dest="normBW", default="none" , type=str)
    optAtac_args.add_argument("-eg", metavar="--effGS"  , help="Effective genome size for Bigwig creation, for when normBW option is enabled.", dest="effGS" , default="none" , type=str)

    # Add the arguments for basic processing groups
    bscProc_args = parser.add_argument_group('Optional Basic Processing Arguments')
    bscProc_args.add_argument("-od" , metavar="--outputDir" , help="* The directory where the output folder will be created.", dest="outputDir", type=str)
    bscProc_args.add_argument("-rg" , metavar="--refGenome" , help="Reference genome fasta file with path.", dest="refGenome", default="none" , type=str)
    bscProc_args.add_argument("-cn" , metavar="--coreNum"   , help="Set core number for Bowtie2 mapping and/or samtools sort.\n- Default: %(default)s" , dest="coreNum", default="4", type=int)
    bscProc_args.add_argument("-sb" , metavar="--sortBams"  , help="If 'yes', sorts the supplied BAM files.\n- Choices: {'yes','no'}\n- Default: %(default)s", dest="sortBams", default="no", type=str, choices=['no','yes'])
    bscProc_args.add_argument("-bl" , metavar="--blackList" , help="Blacklisted regions filtering after peak calling.\n- Choices: {'yes','no'}\n- Default: %(default)s", dest="blackList", default="no", type=str, choices=['no','yes'])
    bscProc_args.add_argument("-bb", metavar="--blBED"      , help="BED file providing the blacklisted regions to be filtered out." , dest="blBED" , default="none" , type=str)
    bscProc_args.add_argument("-hq" , metavar="--hqReads"   , help="High-quality reads (read filtering). \n- Choices: {'no','fastq','bam','both'}\n- Default: %(default)s" , dest="hqReads", default="no", type=str, choices=['no','fastq','bam','both'])
    bscProc_args.add_argument("-qf" , metavar="--qFilter"   , help="Fastq filtering, where qFilter is the minimum quality wanted. E.g.: 30.\n- Default: %(default)s" , dest="qFilter", default="20", type=int)
    bscProc_args.add_argument("-mf" , metavar="--mapqFilter", help="BAM filtering, where mapqFilt is the minimum MAPQ wanted. E.g.: 20.\n- Default: %(default)s", dest="mapqFilter", default="30", type=int)

    # Basic optional group 
    parser.add_argument("-cf" , metavar="--configFile" , help="Config file with the values to the local parameters", dest="config_file", default=None, type=str)


    # http://thomas-cokelaer.info/blog/2014/03/python-argparse-issues-with-the-help-argument-typeerror-o-format-a-number-is-required-not-dict/	

    # Print the help message only if no arguments are supplied	
    if len(sys.argv)==1:	
        parser.print_help()	
        sys.exit(1)	

    # Save the STDOUT output in a log file	
    if parser.parse_args().projDir:	
        logdir="{0}/interimFiles/logs".format(parser.parse_args().projDir)	
        create_dir(logdir)	
        logfile = "{0}/{1}_annotate_variant_wrapper.log".format(logdir, parser.parse_args().sampleName)	
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
    genomeFastaFile   = args.genomeFastaFile	
    genomeSegDupFile  = args.genomeSegDupFile	
    genomeRepMaskFile = args.genomeRepMaskFile	
    genome5bpGCpcFile = args.genome5bpGCpcFile	
    wesProbesFile     = args.wesProbesFile	
    tumorBamFile      = args.tumorBamFile	
    normalBamFile     = args.normalBamFile	
    projDir           = args.projDir
    projName          = args.projName	
    sampleName        = args.sampleName	

    # User variables	
    interimSampleDir  = "{0}/interimFiles/{1}".format(projDir,sampleName)
    interimFileName   = "{0}/interimFiles/{1}/{1}_Mutect2".format(projDir,sampleName)
    inputVariantFile  = "{0}/annInput/{1}_Mutect2.txt".format(projDir,sampleName)
    create_dir(get_file_info(interimFileName)[0])

    main()
