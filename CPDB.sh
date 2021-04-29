#!/bin/bash

for folder in /students/2020-2021/master/CellPhoneDB/Subsets_Smillie/*; do
	cd $folder
	name="$(basename ${folder})"
	echo "working on ${name}"
	cellphonedb method statistical_analysis metadata.csv . --pvalues-result-name ${name}_testPval.txt --threads 60 --counts-data gene_name --project-name $name
	echo "Next"
done
