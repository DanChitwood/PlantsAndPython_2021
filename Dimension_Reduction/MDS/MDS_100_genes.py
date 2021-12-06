#!/usr/bin/env python
# coding: utf-8

# Kara Dobson

# Creating MDS plots with 100 gene subset
# The MDS coordinates were saved to a separate csv (included towards the bottom) so that I didn't have to re-run the MDS each time I wanted to make a plot
# All code to make the MDS itself could be commented out - the csv w/ the coordinates is in the github repo 
# Also note that this MDS was ran on only 100 genes and specified tissue types - this was to get around the clean tissue types not being available at the time of making this
# MDS seems to take a super long time to run for some reason - not sure if its the MDS itself or if my code is funky that makes it take so long

# Import libraries 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import seaborn as sns


test = pd.read_csv('result_table_100.csv', encoding='latin-1')
test.head()

# changing tissue types to the same name - probably a more succinct way to do this, but doesn't really matter since this is just a test
test['Tissue'] = test['Tissue'].replace('leaf' ,'leaves')
test['Tissue'] = test['Tissue'].replace('Leaf' ,'leaves')
test['Tissue'] = test['Tissue'].replace('Leaves' ,'leaves')
test['Tissue'] = test['Tissue'].replace('leave' ,'leaves')
test['Tissue'] = test['Tissue'].replace('leaf tissue' ,'leaves')
test['Tissue'] = test['Tissue'].replace('rosette' ,'rosette leaves')
test['Tissue'] = test['Tissue'].replace('Rosette leaf' ,'rosette leaves')
test['Tissue'] = test['Tissue'].replace('rosettes' ,'rosette leaves')
test['Tissue'] = test['Tissue'].replace('rozettes' ,'rosette leaves')
test['Tissue'] = test['Tissue'].replace('Rosette' ,'rosette leaves')
test['Tissue'] = test['Tissue'].replace('Rosette leaves' ,'rosette leaves')
test['Tissue'] = test['Tissue'].replace('Shoot' ,'shoot')
test['Tissue'] = test['Tissue'].replace('aerial shoots' ,'shoot')
test['Tissue'] = test['Tissue'].replace('whole shoot' ,'shoot')
test['Tissue'] = test['Tissue'].replace('Root' ,'root')
test['Tissue'] = test['Tissue'].replace('roots' ,'root')
test['Tissue'] = test['Tissue'].replace('inflorescences' ,'inflorescence')
test['Tissue'] = test['Tissue'].replace('Inflorescence meristem' ,'inflorescence')
test['Tissue'] = test['Tissue'].replace('flower buds' ,'flower')
test['Tissue'] = test['Tissue'].replace('flower bud' ,'flower')
test['Tissue'] = test['Tissue'].replace('Flower' ,'flower')
test['Tissue'] = test['Tissue'].replace('in vitro cotyledon' ,'cotyledon')
test['Tissue'] = test['Tissue'].replace('root tip' ,'root')

# keep rows with only certain tissue types for now - this will change with the clean dataset
test2 = test[(test == 'leaves').any(axis=1)|(test == 'root').any(axis=1)|(test == 'shoot').any(axis=1)|(test == 'cotyledon').any(axis=1)|(test == 'inflorescence').any(axis=1)|(test == 'flower').any(axis=1)|(test == 'rosette leaves').any(axis=1)|(test == 'hypocotyl').any(axis=1)]
test2['Tissue'].unique()
test2.shape


# making a variable with tissue names to re-merge with MDS dataset later on
sample_names = test2["Tissue"]
sample_names.shape

# taking random subset of the data to speed up processing - don't need this now that the slurm script is running
#test_tissue = test2.sample(frac=0.4, replace=False, random_state=1)
#sample_names = test_tissue["Tissue"]
#test_tissue.shape


# transform data to contain only genes as columns and tissue types as rows
# selecting all gene expression columns and the "Tissue" column
test_tissue = test2[test2.columns[3:103]]
test_tissue = test_tissue.set_index('Tissue')
test_tissue = test_tissue.rename_axis('').rename_axis('Gene', axis='columns')
test_tissue.head()
test_tissue.shape


# Log transform data - do we need to do this? will the data cleaning group do any transformations?
for c in [c for c in test_tissue.columns if np.issubdtype(test_tissue[c].dtype , np.number)]:
    test_tissue[c] += 1
for c in [c for c in test_tissue.columns if np.issubdtype(test_tissue[c].dtype , np.number)]:
    test_tissue[c] = np.log(test_tissue[c])
test_tissue.head()


# apply MDS to get a two-dimensional dataset - euclidean dissimilarity
mds = MDS(n_components=2, random_state=3)
test_mds = mds.fit_transform(test_tissue)


# make dataframe with MDS values
real_mds_df = pd.DataFrame(data = test_mds, columns = ["MDS1", "MDS2"])
real_mds_df.shape


# check that tissue column is good
sample_names.head()
sample_names.shape

# add tissue names back in
real_mds_df.reset_index(drop=True, inplace=True)
sample_names.reset_index(drop=True, inplace=True)
final_real_df = pd.concat([sample_names, real_mds_df], axis = 1)
final_real_df.head()
#final_real_df.to_csv('MDS_coordinates_tissue.csv')

# read in the coordinates csv - downloaded it in the step above so that every time I re-ran this script, I didn't re-run the MDS
#final_real_df  = pd.read_csv('MDS_coordinates.csv', encoding='latin-1')
final_real_df2 = final_real_df.copy()
final_real_df3 = final_real_df.copy()

# plot
sns.lmplot('MDS1', 'MDS2', height=10, data=final_real_df, hue="Tissue", fit_reg=False)
plt.savefig('100genes_tissue.png', facecolor='white', transparent=False)

# change tissue groupings to above/below
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('leaves' ,'Above')
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('rosette leaves' ,'Above')
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('shoot' ,'Above')
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('root' ,'Below')
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('inflorescence' ,'Above')
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('flower' ,'Above')
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('cotyledon' ,'Above')
final_real_df2['Tissue'] = final_real_df2['Tissue'].replace('hypocotyl' ,'Above')

# plot
sns.lmplot('MDS1', 'MDS2', height=10, data=final_real_df2, hue="Tissue", fit_reg=False)
plt.savefig('100genes_above.png', facecolor='white', transparent=False)

# change tissue groupings to veg/rep
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('leaves' ,'Vegetative')
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('rosette leaves' ,'Vegetative')
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('shoot' ,'Vegetative')
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('root' ,'Root')
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('inflorescence' ,'Reproductive')
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('flower' ,'Reproductive')
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('cotyledon' ,'Vegetative')
final_real_df3['Tissue'] = final_real_df3['Tissue'].replace('hypocotyl' ,'Hypocotyl')

# plot
sns.lmplot('MDS1', 'MDS2', height=10, data=final_real_df3, hue="Tissue", fit_reg=False)
plt.savefig('100genes_veg.png', facecolor='white', transparent=False)
