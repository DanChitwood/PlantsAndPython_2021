#!/usr/bin/env python
# coding: utf-8

## Script for hpcc 

import imp
from openTSNE import TSNE
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data = pd.read_csv("/mnt/home/f0103228/Project_docs/data/num_filt_data_cat.csv") # data_genes is a file with only gene expression information (no metadata) 
metadata = pd.read_csv("/mnt/home/f0103228/Project_docs/data/meta_filt_data_cat.csv") # metadata is a file which only contains metadata 

x = data.to_numpy()

# First parameters for tsne need to be established. 

tsne30 = TSNE(
    perplexity=30,
    metric="euclidean",
    n_jobs=8,
    random_state=42,
    verbose=True,
)

# Now we apply tsne to the gene expression data

res = tsne30.fit(x)

# Save object in a file for other teams to use. The result from tsne is a numpy array with 2 columns and as many rows as samples there are.  

df_res = pd.DataFrame(res,columns = ['x_axis','y_axis'])
df_res.to_csv('~/Project_docs/results/tsne_df.csv',index = False)

# We select the relevant columns for applying color in the graph

tissue = metadata["Tissue"].to_numpy()
above_below = metadata["AboveBelow"].to_numpy()
veg_rep = metadata["VegetativeRepro"].to_numpy()


# Plot by tissue

plt.figure(figsize = (20,20),dpi = 300)
plot2 = sns.scatterplot(x = res[:,0], y = res[:, 1], hue = tissue, s = 100, palette = "tab20")
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('/mnt/home/f0103228/Project_docs/results/Tissue_allsamples.png')

# Plot by Above/Below

plt.figure(figsize = (20,20),dpi = 300)
plot2 = sns.scatterplot(x = res[:,0], y = res[:, 1], hue = above_below, s = 100, palette = "tab20")
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('/mnt/home/f0103228/Project_docs/results/AboveBelow_allsamples.png')

# Plot by VegetativeReproductive

plt.figure(figsize = (20,20),dpi = 300)
plot2 = sns.scatterplot(x = res[:,0], y = res[:, 1], hue = veg_rep, s = 100, palette = "tab20")
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('/mnt/home/f0103228/Project_docs/results/VegRepro_allsamples.png')

