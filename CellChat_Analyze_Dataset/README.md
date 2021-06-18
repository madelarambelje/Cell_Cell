# Analyzing a single dataset using Cellchat
### Scripts are based on the original CellChat tutorial: 

- Run "Analyze_Seurat_Object_with_Cellchat.R"
- Select which dataset you want to analyze in line 12 (Info about the datasets are in the main branch)
- Select whether to use permutated dataset if dataset == CellChat in line 15

#### Visualizing:
- If the dataset is the CellChat Tutorial dataset, run visualizations in this script. 


##### - If the dataset is the Lupus Dataset (used in an article that uses ICELLNET): 
- Run "ICELLNET_lupusdata_visualize_cellchat.R"
##### - If the dataset is the BLCA Dataset (used in an article that uses CellPhoneDB): 
- Run "CellPhoneDB_chenzdata_visualize_cellchat.R"


See "CreateSeuratObjectTutorialData.R" to see how the data(sub)set was created from the original data provided in the CellChat tutorial (https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html)
