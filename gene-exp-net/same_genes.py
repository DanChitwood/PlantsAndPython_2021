### our filtered subsets have different genes kept (after the low expression cut) so we need to select only genes that are present in both filtered subsets
import pandas as pd
subset_a=pd.read_csv('rootsub_filtered.csv', header=0, index_col=0)
subset_b=pd.read_csv('leafsub_filtered.csv', header=0, index_col=0)

print("subset a shape is", subset_a.shape) #what do these look like to start
print("subset b shape is", subset_b.shape)

samegenes = subset_a.index.intersection(subset_b.index) # make variable of conserved gened between the two filtered subsets
subset_a = subset_a[subset_a.index.isin(samegenes)] #reduce to only those genes
subset_b = subset_b[subset_b.index.isin(samegenes)]
print(len(samegenes), "genes conserved,\nsubset a new shape is",subset_a.shape, "subset b new shape is", subset_b.shape) # make sure thins look okay

subset_a.to_csv('root_consensus.csv') # make the csv
subset_b.to_csv('leaf_consensus.csv')
