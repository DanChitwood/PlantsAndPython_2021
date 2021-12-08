#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

file = pd.read_csv("tissue_type_dataframe_v2.csv")
catfilter= [col for col in file if not col.startswith('AT')]
file[catfilter]=file[catfilter].replace(np.nan,"Other_NA")
file.to_csv("tissue_type_df_v2_noblank.csv")