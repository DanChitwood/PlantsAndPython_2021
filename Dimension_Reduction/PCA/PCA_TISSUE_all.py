import numpy as np
import pandas as pd

counts = pd.read_csv('/mnt/home/f0103229/Test_project/num_filt_data_cat.csv')

counts_t = counts.transpose()

for c  in counts_t.columns: 
    counts_t[c] += 1
    
for c  in counts_t.columns: 
    counts_t[c] = np.log(counts_t[c])
    
counts_t2 = counts_t.transpose()

from sklearn.decomposition import PCA

pca = PCA(n_components=2)

real_PCs = pca.fit_transform(counts_t2)
real_PCs_df = pd.DataFrame(data = real_PCs, columns = ['PC1', 'PC2'])

sample_annot = pd.read_csv('/mnt/home/f0103229/Test_project/meta_filt_data_cat.csv')

final_real_df = pd.concat([sample_annot, real_PCs_df], axis =1)


import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(15,15))
sns.scatterplot(x='PC1', y='PC2', data=final_real_df, hue='Tissue', legend="full")
plt.legend(loc='center left', bbox_to_anchor=(0.8, 0.5))
plt.savefig('PCA_TISSUE_all.png')


