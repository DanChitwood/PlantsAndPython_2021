#!/usr/bin/env python
# coding: utf-8

# Kara Dobson
# Note: I never got this script to run properly
# Code here is almost the same as the MDS_100_gene.py script, but some things changed for the full dataset

# Import libraries 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import seaborn as sns

# Import dataset
test = pd.read_csv('/mnt/ufs18/rs-008/HRT841_F21/Class_project/Arabidopsis_expression_matrix/DataClean/tissue_type_dataframe_v1.csv', encoding='latin-1')


# Keep rows with only certain tissue types for now - this will change with the clean dataset
# Commented out the code below because the data cleaning group has a complete dataset now (which is the csv loaded in above)
#test2 = test[(test == 'leaves').any(axis=1)|(test == 'root').any(axis=1)|(test == 'whole plant').any(axis=1)|(test == 'seedlings').any(axis=1)|(test == 'rosette leaf').any(axis=1)]
#test2['Tissue'].unique()
#test2.shape


# Making a variable with tissue names to re-merge with MDS dataset later on - never got around to testing this w/ veg/rep or above/below
sample_names = test["Tissue_type"]
sample_names_veg = test["VegetativeRepro"]
sample_names_above = test["AboveBelow"]
sample_names.shape


# Taking random subset of the data - can comment this out once slurm is running
# commented out now because I got slurm to run, but could use these random subsets if its taking too long to run
#test_tissue = test2.sample(frac=0.4, replace=False, random_state=1)
#sample_names = test_tissue["Tissue"]
#test_tissue.shape


# Transform data to contain only genes as columns and tissue types as rows
# First, a dataframe with only gene expression columns & remove descriptor columns
test_tissue = test[test.columns[11:37346,37348]] #first range: genes, second range: tissue type column
# the line above fails; I was trying to select all gene columns and one other column outside that range (the tissue column) but didn't figure it out
test_tissue = test_tissue.set_index('Tissue_type')
test_tissue = test_tissue.rename_axis('').rename_axis('Gene', axis='columns')
test_tissue.head()
#test_tissue.shape


# Log transform data - do we need to do this? will the data cleaning group do any transformations? kept this in for now
for c in [c for c in test_tissue.columns if np.issubdtype(test_tissue[c].dtype , np.number)]:
    test_tissue[c] += 1
for c in [c for c in test_tissue.columns if np.issubdtype(test_tissue[c].dtype , np.number)]:
    test_tissue[c] = np.log(test_tissue[c])
test_tissue.head()


# Apply MDS to get a two-dimensional dataset (euclidean dissimilarity) - this is the step that clogs up memory & takes a long time
mds = MDS(n_components=2, random_state=3)
test_mds = mds.fit_transform(test_tissue)


# Make dataframe with MDS values
real_mds_df = pd.DataFrame(data = test_mds, columns = ["MDS1", "MDS2"])
real_mds_df.shape
real_mds_df.to_csv('MDS_coordinates_full.csv', encoding='utf-8') # saving MDS coordinates so that I don't need to re-run this each time


# Check that tissue column is good
sample_names.head()
sample_names.shape


# Add tissue names back in
real_mds_df.reset_index(drop=True, inplace=True)
sample_names.reset_index(drop=True, inplace=True)
final_real_df = pd.concat([sample_names, real_mds_df], axis = 1)
final_real_df.head()


# Plot
sns.lmplot('MDS1', 'MDS2', height=10, data=final_real_df, hue="Tissue_type", fit_reg=False)
plt.savefig('MDS_full.png', facecolor='white', transparent=False)

