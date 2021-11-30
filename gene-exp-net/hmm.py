## was that ~values as dates~ issue in the original data?

import pandas as pd

pd.set_option('display.max_rows', None)
data=pd.read_csv('gene_FPKM_200501.csv')
print(data.dtypes)
