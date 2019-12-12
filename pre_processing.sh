#!/bin/bash

in=$1

genbank_ncbi(){

	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
	sed 's/ /_/g' assembly_summary_genbank.txt -i

}

idtaxo(){

	mkdir -p tmp
	## remove old files
	dag=$(date | awk -v OFS="_" '{print $3,$2,$6}')
	find tmp/ -name "latest_taxo" -mtime +30 -type f -execdir rm -- '{}' \;
	find tmp/ -name "lineages*" -mtime +30 -type f -execdir rm -- '{}' \;
	## check taxonomy creater
	if [[ ! -s latest_taxo ]]; then
		if [[ ! -d ncbitax2lin ]]; then
			git clone https://github.com/zyxue/ncbitax2lin.git
		fi
		if [[ ! -s tmp/lineages_$dag ]]; then
			make -C ncbitax2lin
			gunzip -c $(find ncbitax2lin -name "line*gz") > tmp/lineages_$dag
		fi
		sed 's/ /_/g;s/^/|/g;s/,/| /1' tmp/lineages_$dag > latest_taxo
		rm ncbitax2lin/*gz* 
	fi

}

db_filter(){
	#mkdir -p tmp
	## Filter out specific terms 
	sed 's/ /_/g' $in | awk -F '","' '$15=="WGS" || $15=="WGA" {print "|"$3"|",$0}' | awk -F '","' '$17!~"cDNA" && $17!~"MDA" && $17!~"Reduced_Representation" && $10>1000' > tmp/1_method_filtered.csv

	awk '{print $1}' tmp/1_method_filtered.csv | awk '!seen[$0]++' > tmp/all_species_found	

	## Get species name, taxid, loc and link
	awk -F'\t' '{print "|"$8"|",$6"\",\""$17"\",\""$20}' assembly_summary_genbank.txt | awk '!seen[$0]++' > tmp/species_taxid
	
	## Join species name with taxid (exact match)
	LC_ALL=C fgrep -wf tmp/all_species_found tmp/species_taxid | sed 's/^|//g;s/| /\",\"/g' | awk -F'","' '{print "|"$2"|",$0}' > tmp/species_found
	rm tmp/species_taxid

	## Filtering unwanted
	awk '{print $1}' tmp/species_found | awk '!seen[$0]++' > tmp/species_taxid_found

	## Link taxid to full taxonomy
##F## Future makes targetted search file.
	LC_ALL=C fgrep -f tmp/species_taxid_found latest_taxo | grep "Arthropoda" | sed 's/ /\t/g;s/,/ /g' > tmp/taxo_search

	## Remove based on search terms in remove_list
	LC_ALL=C fgrep -vf remove_list tmp/taxo_search | sed 's/ /,/g' | awk -F, '$2!=""' > tmp/fin_taxo

	## Taxo and species info
	sed 's/ /_/g;s/\",\"/ /3' tmp/species_found | sort | awk '{if(a!=$1) {a=$1; printf "\n%s%s",$0,FS} else {a=$1;$1="@";printf $0 }} END {printf "\n" }' | sed 's/ @ /@ /g;s/@ /@/g;s/ $//g;s/ /\",\"/g;s/|_/| /g' > tmp/colapsed_species

	awk 'NR==FNR{a[$1]=$2;next} $1 in a{print $2"\",\""a[$1]}' tmp/fin_taxo tmp/colapsed_species | awk -F'","' '{print $1,$3,$2,$5,$4}' > tmp/full_info
	
	awk '{print "|"$1"|"}' tmp/full_info > tmp/final_species

	## Overall table
	LC_ALL=C fgrep -wf tmp/final_species tmp/1_method_filtered.csv | awk '{print $2}' > tmp/2_contam_moved.csv 

	echo -n > tmp/3_final.csv
	awk -F'","' '{print $3,$5}' tmp/2_contam_moved.csv | awk '!seen[$0]++' | sort |\
	while read -r species loc; do
		if [[ ! -z ${loc} ]]; then
			awk -F'","' -v var="${species}" -v var2="${loc}" '$3==var && $5==var2' tmp/2_contam_moved.csv | \
			for i in $(cat); do
				## need " at the beginning and end for the correct format
				lnks=$(awk -v var="${species}" -v var2="${loc}" -v OFS="\",\"" '$1==var && $2==var2 {print ",\""$3,$4,$5"\""}' tmp/full_info )
		 		echo $i$lnks >> tmp/3_final.csv
			done
		else
			awk -F'","' -v var="${species}" -v var2="${loc}" '$3==var && $5==""' tmp/2_contam_moved.csv >> tmp/3_final.csv
		fi 
	done
	sed -i -e 's/\r//g' tmp/3_final.csv

}

idtaxo
db_filter
