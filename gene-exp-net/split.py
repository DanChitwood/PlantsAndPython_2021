### this needs work.. just trying to split the extended metadata off for ease of subsetting & other non-count-based transformations

import pandas as pd

df = pd.read_csv('tissue_type_dataframe_v1.csv', low_memory=FALSE) #load the conjoined meta/data
print(df.shape) #print shape lazily
meta = df['Sample' , 'Project', 'SampleName', 'PMID', 'Genotype', 'Ecotype', 'Tissue', 'TotalReads', 'UniqueMappedRate', 'ReleaseDate', 'Tissue_type', 'VegetativeRepro', 'AboveBelow'] #pull out the meta
print(meta.shape) #lazy

