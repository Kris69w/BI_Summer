ls -t /data/SKDngs01/NGS/OmicsoftRuns/MS00011.PBMC.Guarnieri.sc/Human.B37.3_OmicsoftGene20130723/*.bam |\
while read i; 
do
	name=$(basename "$i" .bam) ;
	dir="/data/SKDngs01/pg/collaborators/junke/SingleCell/PBMC_treat+ctrl" ;
	file="${dir}/bin/sh_${name}.sh"
	echo "generating script for ${name}, basedir ${dir}"  ;
/bin/cat <<EOM > ${file}
#!/bin/bash
#$ -l mem_free=8G
#$ -o ${dir}/log
#$ -e ${dir}/log
#$ -wd ${dir}
### load env from .bashrc
source ~/.bashrc

geneBody_coverage.py -r /data/SKDngs01/pg/annotation/HouseKeepingGenes.bed/hg19.HouseKeepingGenes.nochr.bed -i $i -o ${dir}/plot/o_$name 

EOM
	echo -e "\tqsub file: ${file}\n"  ;
	qsub ${file}
done


