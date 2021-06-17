# Cell-Cell Communication Analysis using scRNA-seq Data
This repository contains the comparison of three novel cell-cell interaction tools. Here we will compare CellPhoneDB, CellChat and ICELLNET. 

These tools predict interactions between cells using ligand and receptors. Each tool in their own way. CellPhoneDB and CellChat are performing permutation testing. ICELLNET calculates the intensity of the ligand-receptor pair and normalizes all the interactions found.


## scRNA-seq Data used

##### CellPhoneDB
Article: [Article data used from Chen et al, Nature Comminications 2020](https://europepmc.org/article/pmc/pmc7545162#MOESM1)

Actual link to data: [CellPhoneDB Article Data Used](https://europepmc.org/articles/PMC7545162/bin/41467_2020_18916_MOESM2_ESM.zip)
#### CellChat
Article: [Tutorial data retrieved from Jin et al., Nature Communications, 2021](https://doi.org/10.1038/s41467-021-21246-9) 

Actual link to data: [CellChat Tutorial Data Used](https://drive.google.com/file/d/1TYsGayCEFa1E2B6HbBb5ctoqMn_ejp4z/view?usp=sharing)
#### ICELLNET
Article: [Tutorial data retrieved from Arazi, et., Nature Immunology 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6726437/)

To acces this data, permission is required from the author.


### Ulcerative Colitus
In this study data additional data is used from [Smillie et al., CellPress 2019](https://pubmed.ncbi.nlm.nih.gov/31348891/), however this data is not further used in the article written.

Actual link to Data: [Smillie Data](https://www.dropbox.com/sh/dn4gwdww8pmfebf/AACXYu8rda5LoLwuCZ8aZXfma?dl=0)

### Comparing 
These tools were tested on their capability in showing relevant results using their visualization options.
Tutorial or article data provided by these tools were used interchangeably on each of the tools, to validate if the tools are able to replicate their results. 
Permutation was performed to detect false postives. Here, we permuted gene, cell labels. Furthermore, we also combined the two permutations, to be sure our data is not making sense, hence no interactions were not expected to be detected. Furthermore, we evalutated the capabilities on showing relevant information using their visualization options of each tool. Also user convinience was tested, by checking how user-friendly the respective tool was. Vignettes were produced for future users.



