Author: Shreya Dey

------
# Cell-Cell communication analysis with ICELLNET
------

This repository explains the use of ICELLNET tool in the research project: "Cell-Cell communication analysis using single-cell RNA seq data: A comparison of three novel tools"

## Contents

- [Introduction](#introduction) 
- [Process](#workflow)
- [Installation of ICELLNET package](#install_ICELLNET)
- [Visualizations](#visualization)
- [Case Study](#casestudy)



## Introduction <a name="introduction"></a>
**ICELLNET** is a transcriptomic-based framework recently developed by [soumelis lab](https://github.com/soumelis-lab) for the analysis of cell to cell communication based on expression of ligand-receptor pairs. In this project ICELLNET was selected as one of the three tools to find a gold standard in analyzing single-cell RNA-seq data. ICELLNET package allows to calculate interaction score based on the interaction pairs present in ligands of central cell to its receptors of cognate peripheral cells. On the basis of calculated interaction scores different visualizations are obtained. Along with the Lupus dataset which is used in tutorial 2 of ICELLNET, ICELLNET was applied on other 2 datasets which were used in the other two tools CellPhoneDB and CellChat.


## Process <a name="workflow"></a>
In this section we have defined details about implementation of ICELLNET on the datasets.

1. Selecting the datasets for comparing the results. 

2. Apply ICELLNET pipeline on the datasets. <a name="score"></a>
Compute interaction score matrix for both inward and outward direction by selecting a central cell from the dataset. We used icellnet.score() function to calculate the interaction score matrix by predefining the direction of communication once as "out" for outward communication and "in" for inward communication.
Only one cell type should be used as a central cell for calculation of interaction score at a time. However, the Partner cell can be a list of multiple cell types. 

3. Use the interaction score matrix to generate the visualization. 


## Installation of ICELLNET package <a name="install_ICELLNET"></a>
The necessary prerequisites to use ICELLNET can be downloaded from this link [prerequisite.r](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Case%20Study/prerequisite.r)


## Visualizations  <a name="visualization"></a>
ICELLNET provides different visualizations in form of network plot, bar plot, bubble plot and heatmap.
Predefined functions were used to generate the plots.

1. Network Plot

Function network.create() was used to create [network plots](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Visualizations/networkplot.png) to show communication between a central cell to a selection of multiple peripheral cells. The interaction score is normalized on a scale of 1-10 before plotting it into a network graph.Edges represent cell types and nodes represent communication between the cell types. The arrows show in which direction the communication is present and the width of arrows represents the intensity of communication based on the normalized score (1-10). In the plot, **R** and **L** denotes **receptors** and **ligands** expressed by central cell respectively. LC(R)-in signifies communication is from ligands of peripheral cells to receptors expressed by central cell LC. LC(L)-out signifies communication is from ligands of central cell LC to receptors of peripheral cells.


2. Bar Plot

To calculate the contribution of each family of molecules to the global communication scores, function LR.family.score() was used. The ligand-receptor interaction score calculated previously in [Process step 2](#score)  is passed as an argument. When the plot argument is set as FALSE it creates a table and displays it in the R console and when it is set to TRUE it creates a [barplot](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Visualizations/barplot.png). This plot helps in understanding the distribution of communication molecules based on families or subfamilies of molecules which are present in the ICELLNET database. The category defined as other is not yet classified in the ICELLNET database. Hence, it is difficult to say which families or subfamilies the cells are contributing.


3. Bubble Plot/Bubble Plot

Bubble Plot in ICELLNET shows a score on interaction between the most contributing ligand/receptor pairs to the communication score. This allows to identify specific individual interactions that can drive the intercellular communication and should be confirmed experimentally. LR.balloon.plot() function was used to create [bubble plots](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Visualizations/balloon_plot.png). The contribution scores can be sorted based on 2 criterias by 'sum' or by 'var'. 'sum' is used to highlight highest contribution factors while 'var' is used to highlight the differences. We used sum for our visualizations as we wanted to check which interaction pairs contribute the highest. The color shows family of molecule and the size of bubble represents interaction score on a scale of 0-100. In the plot, only top 30 interactions are shown for cell type LC. Specific interaction pairs can be selected by using LR score matrix which is shown in casestudy 2 and casestudy 3.


4. Heatmap

LR.heatmap() function was used to create [heatmaps](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Visualizations/heatmap.png) based on the calculated interaction score. It shows identical details as the bubble plot but in a heatmap. This function is not yet published by the original author of ICELLNET therefore it is not possible to mention further details of how the function works. As part of our project the source code was exclusively shared with us by the author to create visualizations. Hence,in this project only the final results are published.  


## Casestudy  <a name="casestudy"></a>

This section contains the link for the code of ICELLNET implementation for tutorial datasets.

1) [CellChat tutorial dataset](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Case%20Study/CellChat_Permutations_ICELLNET.R)

2) [CellPhoneDB tutorial dataset](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Case%20Study/Cellphonedb_ICELLNET.R)

3) [ICELLNET tutorial dataset](https://github.com/madelarambelje/Cell_Cell/blob/shreya/Case%20Study/Lupus_ICELLNET.R)

