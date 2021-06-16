#!/bin/bash

# Making Directory for input
mkdir $1

# Creating input from seurat object provided by the user
Rscript rdsToMBF.R $1 $2 

# Move to input directory 
cd $1

# Running statistical method of cellphonedb using the input made previously
cellphonedb method statistical_analysis metadata.csv . --pvalue 0.05 --pvalues-result-name ${1}_testPval.txt --threads $3 --counts-data gene_name

# Moving to output created by the statistical method
cd out/

# Creating heatmap by cellphoneDB and network file for manual visualizations
cellphonedb plot heatmap_plot ../metadata.csv --pvalue 0.05 --pvalues-path ${1}_testPval.txt --count-network-name count_network.txt


