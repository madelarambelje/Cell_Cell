#!/bin/bash

mkdir $1

Rscript rdsToMBF.R $1 $2 

cd $1

cellphonedb method statistical_analysis metadata.csv . --pvalue 0.05 --pvalues-result-name ${1}_testPval.txt --threads $3 --counts-data gene_name

cd out/

cellphonedb plot heatmap_plot ../metadata.csv --pvalue 0.05 --pvalues-path ${1}_testPval.txt --count-network-name count_network.txt

cd out/

Rscript ../../../network_heatmap.R 
