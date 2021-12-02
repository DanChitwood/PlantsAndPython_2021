#!/usr/bin/env python
# coding: utf-8

# ### _Arabidopsis_ gene sequence expression mapper graph with centroid lens
# 
# **Imort useful packages / modules**

# In[1]:


# import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ML tools
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import MinMaxScaler
from sklearn import ensemble

# For output display
from IPython.display import IFrame

# # If running locally, set current directory as projdir
# projdir = '.'


# **For google collab**

# In[2]:


# # Only if running in Google Colab..!!
# # DO NOT run this cell if running locally - simply comment it out.
# from google.colab import drive
# drive.mount('/content/gdrive/')

# projdir = '/content/gdrive/MyDrive/PlantsAndPython-2021-10-22'
# sys.path.append(projdir)


# **Modules from Sourabh Palande files**

# In[2]:


# # import helper_functions
# from helper_functions import loaddata # for data loading
from helper_functions import colorscale_from_matplotlib_cmap # for kmapper color palette

# # import Nicolaou et al. 2011 lense function
# from lenses import fsga_transform

# keppler mapper
import kmapper as km


# **Import database**

# In[3]:


df = pd.read_csv("37336genes_11317samples_20tissues_tSNE-centroild-lens.csv")
# This database has:
# - tissue types labeled according to the hypothesis generation group classification
# - euclidean distances of each coordinate t-SNE point to its respective tissue type centroid
# - 37,336 genes from the tissue_type_dataframe_v1.csv file in the CleanData HPCC directory
# - 11,316 Arabidopsis samples from the tissue_type_dataframe_v1.csv file in the CleanData HPCC directory
# This database was made 01/Dec/21

# NOTE: It takes about 30 minutes to load using 1 core, 1 node and 16GB in the HPCC


# **Exploring database**

# In[4]:


df_name = 'Centroid_Lens_Database'
print("The current dataframe name is:", df_name)
print("rows, columns =", df.shape)
print("number of elements =", df.size)


# **Subsetting gene expression columns**

# In[79]:


genes = list(df.columns[11:-30]) # create list with the gene names (first gene is column 11 and last gene is column -30)
len(genes) # check how many genes you're using (maximum is 37,336)


# **Set factors and factors levels**

# In[73]:


factors = ['Tissue','VegetativeRepro','AboveBelow','Sample Type']
levels = ['Root','Root','Below','knl2 mutant line (flowering buds)']

filter_by_factor, filter_by_level = ('Tissue', 'Root')
# filter_by_factor, filter_by_level = ('VegetativeRepro', 'Root')
# filter_by_factor, filter_by_level = ('AboveBelow', 'Below')
# filter_by_factor, filter_by_level = ('SampleName', 'knl2 mutant line (flowering buds)')

# color_by_factor, color_by_level = ('Tissue', 'Root')
# color_by_factor, color_by_level = ('VegetativeRepro', 'Root')
color_by_factor, color_by_level = ('AboveBelow', 'Below')


# **Initialize a KeplerMapper object**

# In[74]:


# Initialize mapper object
mymapper = km.KeplerMapper(verbose=1)

# Define Nerve
nerve = km.GraphNerve(min_intersection=1)


# **Define lens**
# 
# According to Dan's description: _"take the centroid/median of each tissue cluster, and the lens is the eucledian distance of each sample to its respective tissue center"_

# In[75]:


# Centroid lens
Clens = df["eucl_dist"] # the euclidean distances are found in the "eucl_dist" column
lens_type = 'Centroid'
#plt.plot(Clens) # plot the lens to see how well represents the data


# **Define cover:**
# 
# Overlap must be between 0 and 100. Intervals must be less than 90: try between 25 to 85.

# In[76]:


# Define cover
cubes, overlap = (100, 75) # cubes = intervals
cover = km.cover.Cover(n_cubes=cubes, perc_overlap=overlap/100.)


# **Define clustering algorithm:**
# 
# DBSCAN with default parameters. Metric: correlation distance (1 - correlation) between a pair of gene expression profiles.

# In[77]:


# Define clustering algorithm
clust_metric = 'correlation'
clusterer = DBSCAN(metric=clust_metric)


# **Construct the mapper graph:**
# 
# Keep an eye on the number of hypercubes, nodes and edges reported by the algorithm. You can change the graph size by changing the cover parameters.

# In[78]:


# Create mapper 'graph' with nodes, edges and meta-information.
graph = mymapper.map(lens=Clens,
                     X=df[genes],
                     clusterer=clusterer,
                     cover=cover,
                     nerve=nerve,
                     precomputed=False,
                     remove_duplicate_nodes=True)


# **Kmapper coloring**

# In[80]:


# Color nodes by specified color_by_factor, color_by_level

df[color_by_factor] = df[color_by_factor].astype('category')
color_vec = np.asarray([0 if(val == color_by_level) else 1 for val in df[color_by_factor]])
cscale = colorscale_from_matplotlib_cmap(plt.get_cmap('coolwarm'))


# **Set coloring levels as kmapper tooltips**

# In[81]:


# show color_by_factor levels in tooltip

temp = ['({}, {})'.format(str(p[0]), str(p[1])) for p in zip(df[color_by_factor], df[filter_by_factor])]
df['tooltips'] = temp


# **Create and save kmapper graph as html**

# In[ ]:


# Specify file to save html output
fname = 'LensType_{}_ColorBy_{}_Tips_{}_Data_{}_Cubes_{}_Overlap_{}_Genes_{}.html'.format(lens_type,
                                                              color_by_factor,
                                                              filter_by_factor,
                                                              df_name,
                                                              cubes,
                                                              overlap,
                                                              len(genes))

figtitle = 'Lens type: {}, Tips {}, Color by {} ({}), Database: {}, intervals {}, overlap {}, genes {}'.format(lens_type,
                                                                                                              filter_by_factor,
#                                                                                                               filter_by_level,
                                                                                                              color_by_factor,
                                                                                                              color_by_level,
                                                                                                              df_name,
                                                                                                              cubes, 
                                                                                                              overlap/100.0,
                                                                                                              len(genes))

fpath = '/mnt/home/f0103237/SLURM_mapper_outputs/' + fname # is this synthax correct if I run it in the HPCC?

# Create visualization and save to specified file
_ = mymapper.visualize(graph,
                       path_html=fpath,
                       title=figtitle,
                       color_values=color_vec,
                       color_function_name=color_by_factor,
                       colorscale=cscale,
                       custom_tooltips=df['tooltips'])

# Load the html output file
IFrame(src=fpath, width=1000, height=800)

