## was that ~values as dates~ issue in the original data?
## i don't think so
## yes i use git improperly, but i find the web gui nice for writing code in ¯\_(ツ)_/¯

import pandas as pd

pd.set_option('display.max_rows', None)
data=pd.read_csv('gene_FPKM_200501.csv')
print(data.dtypes)
