#!/usr/bin/env python
# coding: utf-8

# ### _Arabidopsis_ gene sequence expression mapper graph with centroid lens

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


# In[3]:


# # Make sure your last path is the one were you have this script and your data
# sys.path


# **Modules from Sourabh's files**

# In[4]:


# # import helper_functions
# from helper_functions import loaddata # for data loading
from helper_functions import colorscale_from_matplotlib_cmap # for kmapper color palette

# # import Nicolaou et al. 2011 lense function
# from lenses import fsga_transform

# keppler mapper
import kmapper as km


# **Import database**

# In[5]:


centroids = pd.read_csv("./metadata_centroids_1000genes_19642samples.csv")
# This database contains 1000 genes from the MapRateFiltered_v1.csv file in the CleanData HPCC directory
# This database has tissue types labeled according to the hypothesis generation group classification
# This database was made 24/Nov/21


# **Defining database as df**

# In[6]:


df = centroids # set dataframe
df_name = 'centroids'
print("The current dataframe name is:", df_name)
print("rows, columns =", df.shape)
print("number of elements =", df.size)


# **Subsetting gene expression columns**

# In[25]:


genes = list(df.columns[40:540]) # create list with the gene names
len(genes)


# **Set factors and factors levels**

# In[8]:


factors = ['Tissue','VegetativeRepro','AboveBelow']
levels = ['Root','Root','Below']

# filter_by_factor, filter_by_level = ('Tissue', 'Root')
# filter_by_factor, filter_by_level = ('VegetativeRepro', 'Root')
# filter_by_factor, filter_by_level = ('AboveBelow', 'Below')

color_by_factor, color_by_level = ('Tissue', 'Root')
# color_by_factor, color_by_level = ('VegetativeRepro', 'Root')
# color_by_factor, color_by_level = ('AboveBelow', 'Below')


# **Initialize a KeplerMapper object**

# In[9]:


# Initialize mapper object
mymapper = km.KeplerMapper(verbose=1)

# Define Nerve
nerve = km.GraphNerve(min_intersection=1)


# **Define lens**
# 
# According to Dan's description: _"take the centroid/median of each tissue cluster, and the lens is the eucledian distance of each sample to its respective tissue center"_

# In[10]:


# Centroid lens
Clens = df["eucl_dist"] # the euclidean distances are found in the "eucl_dist" column
lens_type = 'Centroid'
# plt.plot(Clens) # plot the lens to see how well represents the data


# In[ ]:


# # Multiple centroids lens
# CNTRDSlens = df.iloc[:, 20:22]
# lens_type = "MultipleCentroids"
# plt.plot(CNTRDSlens)


# **Define cover:**
# 
# Overlap must be between 0 and 100. Intervals must be less than 130.

# In[18]:


# Define cover
cubes, overlap = (100, 90) # cubes = intervals
cover = km.cover.Cover(n_cubes=cubes, perc_overlap=overlap/100.)


# **Define clustering algorithm:**
# 
# DBSCAN with default parameters. Metric: correlation distance (1 - correlation) between a pair of gene expression profiles.

# In[19]:


# Define clustering algorithm
clust_metric = 'correlation'
clusterer = DBSCAN(metric=clust_metric)


# **Construct the mapper graph:**
# 
# Keep an eye on the number of hypercubes, nodes and edges reported by the algorithm. You can change the graph size by changing the cover parameters.

# In[20]:


# Create mapper 'graph' with nodes, edges and meta-information.
graph = mymapper.map(lens=Clens,
                     X=df[genes],
                     clusterer=clusterer,
                     cover=cover,
                     nerve=nerve,
                     precomputed=False,
                     remove_duplicate_nodes=True)


# **Kmapper coloring**

# In[21]:


# Color nodes by specified color_by_factor, color_by_level

df[color_by_factor] = df[color_by_factor].astype('category')
color_vec = np.asarray([0 if(val == color_by_level) else 1 for val in df[color_by_factor]])
cscale = colorscale_from_matplotlib_cmap(plt.get_cmap('coolwarm'))


# **Set coloring levels as kmapper tooltips**

# In[22]:


# show color_by_factor levels in tooltip

temp = ['({})'.format(str(p[0])) for p in zip(df[color_by_factor])]
df['tooltips'] = temp


# **Create and save kmapper graph as html**

# In[ ]:


# Specify file to save html output
fname = 'LensType_{}_ColorBy_{}_Data_{}_Cubes_{}_Overlap_{}_Genes_{}.html'.format(lens_type,
#                                                               filter_by_factor,
                                                              color_by_factor,
                                                              df_name,
                                                              cubes,
                                                              overlap,
                                                              len(genes))

figtitle = 'Lens type: {}, Color by {} ({}), Database: {}, intervals {}, overlap {}, genes {}'.format(lens_type,
#                                                                                                               filter_by_factor,
#                                                                                                               filter_by_level,
                                                                                                              color_by_factor,
                                                                                                              color_by_level,
                                                                                                              df_name,
                                                                                                              cubes, 
                                                                                                              overlap/100.0,
                                                                                                              len(genes))

fpath = './' + fname # is this synthax correct if I run it in the HPCC?

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


# In[ ]:




