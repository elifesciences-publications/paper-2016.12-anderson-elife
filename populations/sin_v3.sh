#!/bin/bash


### get paths to bam files
	mkdir ~/sin/
	ls */mapping/*_InDel.bam | tr '\t' '\n' > ~/sin/files
	sed  -i 's/^/\/mnt\/icy_4\/droseu\/mapped\//g' ~/sin/files
	
### call indels genome-wide
	### run VarScan for indels
		chrInd () {
		chr=${1}
		samtools mpileup \
		--fasta-ref /mnt/icy_4/droseu/pipeline/reference/Dmel_6.04_hologenome_v2.fasta \
		--region ${chr} \
		--bam-list ~/sin/files | \
		java -jar ~/VarScan.v2.3.9.jar \
		mpileup2indel > ~/sin/pooled/indels_${chr}.vcf

		}
		export -f chrInd

		nohup parallel --gnu -j5 chrInd ::: 2L 2R 3L 3R X & 
	
	### parse vcf files
		parseVCF () {

		chr=${1}

		sampNames=$( cat ~/sin/files | cut -f6 -d'/' | tr '\n' ',' )

		cat ~/sin/pooled/indels_${chr}.vcf | grep -v "Chrom" | grep -m1 "47507" | awk -v samps="${sampNames}" '
		BEGIN {
		print "chr\tpos\tref\talt\tpop\tnAlt\trd"
		nSamps=split(samps, sampSplit, ",")
		}
		{
		for(i=11; i<=NF; i++) {
		printf $1"\t"$2"\t"$3"\t"$4"\t"sampSplit[(i-10)]"\t"
		split($i, sp, ":")
		print sp[4]"\t"sp[2]
		}
		}' #> cat ~/sin/pooled/indels_${chr}.delim

		}
		export -f parseVCF 

		parseVCF 2L > ~/sin/pooled/indels_2L.delim
		parseVCF 2R > ~/sin/pooled/indels_2R.delim
		parseVCF 3L > ~/sin/pooled/indels_3L.delim
		parseVCF 3R > ~/sin/pooled/indels_3R.delim
		parseVCF X > ~/sin/pooled/indels_X.delim

		parallel --gnu -j5 parseVCF ::: 2L 2R 3L 3R X
	
### download DPGP data
	mkdir ~/sin/dpgp
	cd ~/sin/dpgp
	
	wget http://pooldata.genetics.wisc.edu/DSPR_round1_indels.tar.bz2
	wget http://pooldata.genetics.wisc.edu/dpgp2_round1_indels.tar.bz2
	wget http://pooldata.genetics.wisc.edu/dpgp3_round1_indels.tar.bz2
	wget http://pooldata.genetics.wisc.edu/dgrp_round1_indels.tar.bz2
	wget http://pooldata.genetics.wisc.edu/DSPR_round1_indels.tar.bz2
	wget http://pooldata.genetics.wisc.edu/ages_round1_indels.tar.bz2
	wget http://pooldata.genetics.wisc.edu/DGN11/BERGMAN_round1_indels.tar.gz
	wget http://pooldata.genetics.wisc.edu/DGN11/CLARK_round1_indels.tar.gz
	wget http://pooldata.genetics.wisc.edu/DGN11/NUZHDIN_round1_indels.tar.gz
	wget http://pooldata.genetics.wisc.edu/DGN11/POOL_round1_indels.tar.gz
	
	tar xvfj DSPR_round1_indels.tar.bz2
	tar xvfj dgrp_round1_indels.tar.bz2
	tar xvfj dpgp2_round1_indels.tar.bz2
	tar xvfj dpgp3_round1_indels.tar.bz2
	tar xvfz BERGMAN_round1_indels.tar.gz
	tar xvfz CLARK_round1_indels.tar.gz
	tar xvfz NUZHDIN_round1_indels.tar.gz
	tar xvfz POOL_round1_indels.tar.gz

### bgzip and tabix

zip_index () {
echo ${1}
bgzip -@20 -c ${1} > ${1}.gz
tabix -p vcf ${1}.gz
}

export -f zip_index
parallel --gnu -j1 zip_index ::: $( ls /home/bergland/sin/dpgp/NUZHDIN_round1_indels )

parallel --gnu -j1 zip_index ::: $( ls /home/bergland/sin/dpgp/BERGMAN_round1_indels )
parallel --gnu -j1 zip_index ::: $( ls /home/bergland/sin/dpgp/POOL_round1_indels )
parallel --gnu -j1 zip_index ::: $( ls /home/bergland/sin/dpgp/CLARK_round1_indels )


### extract out indel polymorphism data from dpgp data
getIndelCall () {	
tabix -p vcf ${1} 8:12236496-12236496 | nl | awk -v fn=${1} '{
split($9, sp, ";")
split(sp[2], spp, "=")
af=spp[2]


}
END {
fn_l=split(fn, fns, "\\/")

if(NR==0) print fns[fn_l]",0"
if(NR>0) print fns[fn_l]","af
}' | sed 's/_INDELS.vcf.gz//g'
}
export -f getIndelCall

parallel --gnu -j1 getIndelCall ::: $( ls /home/bergland/sin/dpgp/*.vcf.gz ) > /home/bergland/sin/dpgp_freqs.csv
parallel --gnu -j1 getIndelCall ::: $( ls /home/bergland/sin/dpgp/BERGMAN_round1_indels/*.vcf.gz ) > /home/bergland/sin/Bergman_dpgp_freqs.csv
parallel --gnu -j1 getIndelCall ::: $( ls /home/bergland/sin/dpgp/POOL_round1_indels/*.vcf.gz ) > /home/bergland/sin/Pool_dpgp_freqs.csv
parallel --gnu -j1 getIndelCall ::: $( ls /home/bergland/sin/dpgp/CLARK_round1_indels/*.vcf.gz ) > /home/bergland/sin/Clark_dpgp_freqs.csv
	


#### vcf to thap
cat /mnt/icy_2/indoorAse/DGRP_prior/dgrp2.vcf | grep -v "#" | awk '{
for(i=10; i<=NF; i++) {
if($i=="0/0") { 
printf $4" " > "/home/bergland/sin/dgrp."$1".thap"
} else if($i=="1/1") {
printf $5" " > "/home/bergland/sin/dgrp."$1".thap"
} else {
printf "N " > "/home/bergland/sin/dgrp."$1".thap"
}
}
printf "\n" > "/home/bergland/sin/dgrp."$1".thap"
print $3" "$1" "$2" "$4" "$5 > "/home/bergland/sin/dgrp.inp"
}'







