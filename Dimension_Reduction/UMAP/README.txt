This folder contains all files related to UMAP within the Dimension Reduction group.
Author: Robin Waterman

-------------------------------------------------------------------------
Subsets : folder contains all files related to subsets of the data

	UMAP_3gene.ipynb : for 3-gene subset, downloaded directly from the website, called "results_table.csv" in this folder
	
	UMAP_tiny75.ipynb : for 1000-gene subset called "MapRateFiltered_tiny_75.csv" in this folder
	
	UMAP_1000Gene.html, UMAP_1000Gene_2D.png, UMAP_1000Gene_3D.png : figure outputs from UMAP_tiny75 script

------------------------------------------------------------------------
Main: folder contains all files related to the full dataset

	NoNA.py : python script for changing blank cells (NaN when read into pandas) in the categorical columns of "tissue_type_dataframe_v2.csv" into the string "Other_NA" and exporting as "tissue_type_df_v2_noblank.csv". Without this step, plotly fails to color by categories.
	*note: used sed -i 's/Tissue.1/Tissue_type/g' tissue_type_dataframe_v2.csv in bash first because plotly gave an error with Tissue.1

	UMAP_All.ipynb : Main Jupyter notebook from which .py scripts were downloaded and used to run UMAP and produce figures. Comments throughout notebook explain steps. Runs on "tissue_type_df_v2_noblank.csv".

	UMAP_All_SLURM : SLURM submission script for UMAP_All.py (downloaded from UMAP_All.ipynb)

	Figures: folder contains all figures produced from analyses of full dataset
	key to file names:
		TissueType/AboveBelow/VegRepro = color categorization
		10/30 = n_neighbors parameter
		.1/.05 = min_dist parameter
		.png = static (2D or 3D)
		.html = 3D interactive
		v1/v2 = version of dataset file used. The main difference is that in v1, all rows wih NA for any of the 3 tissue categories (TissueType, AboveBelow, VegetativeRepro) were removed, while in v2 these were retained. See DataCleaning group for more info. 