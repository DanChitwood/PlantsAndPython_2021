import pandas as pd
import numpy as np

subset = pd.read_csv('rootsub.csv', header=0, index_col=0)

def ex(series, thresh):
    return (series <= thresh).sum() #defines function to sum number of values below threshold in row
thresh=10 #sets threshold to 10
subset['lowexpr']=subset.apply(
    func=lambda row: ex(row, thresh), axis=1) #applies function to all rows & appends output in new column

sampls=len(roots.columns)-1 #obtain number of samples
subset_expr = subset[subset['lowexpr'] < (0.9*sampls)] #remove genes with low expression
subset_expr=subset_expr.drop('lowexpr', 1) #remove expression column

subset_log = subset_expr.apply(lambda x: np.log2(x+1))
subset_log.to_csv('rootsub_filtered.csv')

