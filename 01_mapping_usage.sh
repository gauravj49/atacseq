
# Run in parallel
cmd="parallel ::: "
for s in scripts/01_map_fastQ_bowtie2/AGRad_ATACseq_MUC001/*.sh
do 
	chmod 775 ${s}
	cmd=$(echo "${cmd} ${s}")
done

echo ${cmd}
# to run
eval ${cmd}
