echo "0.000000000001 0 0.000000000001" > range_list
echo "0.0000000001 0 0.0000000001" >> range_list
echo "0.00000001 0 0.00000001" >> range_list
echo "0.000001 0 0.000001" >> range_list
echo "0.00001 0 0.00001" >> range_list
echo "0.0001 0 0.0001" >> range_list
echo "0.001 0 0.001" >> range_list
echo "0.01 0 0.01" >> range_list
echo "0.1 0 0.1" >> range_list
echo "1.0 0 1.0" >> range_list

for chr in 1 3 5 7 9 11 13 15 17 19 21
do
	plink --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb \
	--chr ${chr}  \
	--out parsed_data/rs_chr${chr} --write-snplist
	cat parsed_data/rs_chr${chr}.snplist >> parsed_data/rs_odd.txt
done

for chr in 2 4 6 8 10 12 14 16 18 20 22
do
	plink --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb \
	--chr ${chr}  \
	--out parsed_data/rs_chr${chr} --write-snplist
	cat parsed_data/rs_chr${chr}.snplist >> parsed_data/rs_even.txt
done


### summ stats from: https://alkesgroup.broadinstitute.org/UKBB

for phen in disease_ALLERGY_ECZEMA_DIAGNOSED blood_EOSINOPHIL_COUNT disease_ASTHMA_DIAGNOSED mental_NEUROTICISM
do

	cut -f1 sumstats/${phen}.sumstats | uniq -d > exclude_${phen}    ## filter multi-allelic

	cat parsed_data/rs_odd.txt exclude_${phen} > exclude_odd_${phen}
	cat parsed_data/rs_even.txt exclude_${phen} > exclude_even_${phen}

	awk '{print $1,$10}' sumstats/${phen}.sumstats  | uniq -u -f 1 > tmp.SNP.pvalue_${phen}

	plink --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb \
	--score sumstats/${phen}.sumstats 1 4 8 header \
	--q-score-range range_list tmp.SNP.pvalue_${phen} \
	--exclude exclude_${phen} \
	--out parsed_data/${phen}

	plink --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb \
	--score sumstats/${phen}.sumstats 1 4 8 header \
	--q-score-range range_list tmp.SNP.pvalue_${phen} \
	--exclude exclude_odd_${phen} \
	--out parsed_data/${phen}_even

	plink --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb \
	--score sumstats/${phen}.sumstats 1 4 8 header \
	--q-score-range range_list tmp.SNP.pvalue_${phen} \
	--exclude exclude_even_${phen} \
	--out parsed_data/${phen}_odd

done
rm range_list exclude_* tmp.SNP.pvalue_*
