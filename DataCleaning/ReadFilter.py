# Py Script to be ran with corresponding slurm script, should take about a half hour to run as a job in dev-node
import pandas as pd

# Read meta data and full gene matrix in
df = pd.read_csv('gene_FPKM_200501.csv')
df_meta = pd.read_table('Arabidopsis_metadata.tsv')

# Take a subset of samples that have poor UniqueMappedRate, you could adjust the filter here
bad_df = df_meta[df_meta['UniqueMappedRate'] < 0.75]
bad_samples = list(bad_df['Sample'])

# The gene expression data has sample as columns and genes as index, need to switch
Gene_list = df['Sample']
print('Shape before transpose', df.shape)
df = df.transpose()
print('Shape after transpose', df.shape)

# And lets make columns genes
df.columns = Gene_list
df.drop(labels = 'Sample', axis = 0, inplace=True)
# After transposing the sample names become the index so we need to make column
df['Sample'] = df.index

# Filter samples and record how many were removed
print('Shape before filter', df.shape)
num1 = df.shape[0]
# Use list of bad samples we created earlier
df = df[~df['Sample'].isin(bad_samples)]
print('Shape after filter', df.shape)
num2 = df.shape[0]
print('Samples removed:', num1 - num2)

# Merge filtered data with matching meta data
df = pd.merge(df_meta, df, on='Sample')

# Print out some useful info
print('End shape', df.shape)
print('MapRate Before', df_meta['UniqueMappedRate'].describe())
print('MapRate After', df['UniqueMappedRate'].describe())

# Index = False to prevent "Unnamed: 0" column from appeared in saved file
df.to_csv('MapRateFiltered.csv', index = False)
