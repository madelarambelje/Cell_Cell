Author: Shreya Dey

# Cell-Cell communication analysis with ICELLNET
This repository explains the use of ICELLNET tool in the research project: "Cell-Cell communication analysis using single-cell RNA seq data:A comparison of three novel tools"

## Contents

-[Introduction](#introduction) 
-[Workflow](#workflow)
-[Installation of ICELLNET package](#install_ICELLNET)
-[Visualizations](#visualization)
-[Case Study](#casestudy)

## Introduction {#introduction}
**ICELLNET** is a transcriptomic-based framework recently developed by [soumelis lab](https://github.com/soumelis-lab) for the analysis of cell to cell communication based on expression of ligand-receptor pairs. In this project ICELLNET was selected as one of the three tools to find a *gold standard* in analyzing single-cell RNA-seq data. ICELLNET package allows to calculate interaction score based on the interaction pairs present in ligands of central cell to its receptors of cognate peripheral cells. On the basis of calculated interaction scores different visualizations are obtained. Along with the Lupus dataset which is used in tutorial 2 of ICELLNET, ICELLNET was applied on other 2 datasets which were used in the other two tools CellPhoneDB and CellChat.


## Workflow {#workflow}
In this section we have defined stepwise details about implementation of ICELLNET on the datasets.

#1. Selecting the tutoprial datasets for the tools.

#2. Applying ICELLNET pipeline for all the datasets. {#score}
Compute interaction score matrix for both inward and outward direction by selecting a central cell from the dataset. We used icellnet.score() function to calculate the interaction score matrix by predefining the direction of communication once as "out" for outward communication and "in" for inward communication.
Only one cell type should be used as a central cell for calculation of interaction score at a time. However, the Partner cell can be a list of multiple cell types. 

#3. Use the interaction score matrix to generate the visualization. 


## Installation of ICELLNET package {#install_ICELLNET}
The necessary prerequisites to use ICELLNET can be downloaded from this link [prerequisite.r](https://github.com/madelarambelje/Cell_Cell/blob/shreya/prerequisite.r)


## Visualizations {#visualization}
ICELLNET provides different visualizations in form of network plot, bar plot, bubble plot and heatmap.
Predefined functions were used to generate the plots.

#1. Network Plot
Function network.create() was used to create Network plots to show communication between a central cell to a selection of multiple peripheral cells. The interaction score is normalized on a scale of 1-10 before plotting it into a network graph.Edges represent cell types and nodes represent communication between the cell types. The arrows show in which direction the communication is present and the width of arrows represents the intensity of communication based on the normalized score (1-10)


#2. Bar Plot
To calculate the contribution of each family of molecules to the global communication scores LR.family.score() was used where ligand-receptor interaction score calculated previously in [step 2 of Workflow](#score) is passed as an argument to compute the contribution of families or subfamilies of molecules to the global communication scores. When the plot argument is set as FALSE it creates a table and displays it in the R console and when it is set to TRUE it creates a barplot. This plot helps in understanding the distribution of communication molecules based on families or subfamilies of molecules which are present in the ICELLNET database.

#3. Bubble Plot
Bubble Plot in ICELLNET shows a score on interaction between the most contributing ligand/receptor pairs to the communication score. This allows to identify specific individual interactions that can drive the intercellular communication and should be confirmed experimentally. LR.balloon.plot() function was used to create bubble plots. The contribution scores can be sorted based on 2 criterias by 'sum' or by 'var'. sum is used to highlight highest contribution factors while var is used to highlight the differences. We used sum for our visualizations as we wanted to check whcih interaction pairs contribute the highest.

#4. Heatmap
LR.heatmap() function was used to create heatmaps based on the calculated interaction score. It shows identical details as the bubble plot but in a heatmap. This function is not yet published by the original author of ICELLNET therefore it is not possible to mention further details of how the function works. As a part of our project the source code was exclusively shared with us by the author to create visualizations.Hence,in this vignette only the final results are published.  

## Casestudy {#casestudy}

#1) [ICELLNET tutorial dataset](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Case%20Study/Lupus_ICELLNET.R)

#2) [CellPhoneDB tutorial dataset](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Case%20Study/Cellphonedb_ICELLNET.R)

#3) [CellChat tutorial dataset](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Case%20Study/CellChat_Permutations_ICELLNET.R)

