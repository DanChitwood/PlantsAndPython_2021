#!/usr/bin/env python
# coding: utf-8

# In[5]:


# Step 1: Import packages, functions, and classes
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
import pandas as pd
from matplotlib import pyplot as plt

# Step 2: Get data
dat = pd.read_csv("multipleGenes.csv", index_col=0)
print(dat.head())

xdata = dat[["AT1G22630","AT1G22620","AT1G22610"]]
xdata.head()

y_tissues = dat["Tissue"]
y_stress = dat["Treatment"]
y_genotype = dat['Genotype']

x = np.array(xdata)
y = np.array(y_tissues)
y2 = np.array(y_stress)
y3 = np.array(y_genotype)

# Step 3: Create a model and train it
model = LogisticRegression(solver='liblinear', C=10.0, random_state=0)
model.fit(x, y)

# Step 4: Evaluate the model
p_pred = model.predict_proba(x)
y_pred = model.predict(x)
score_ = model.score(x, y)
conf_m = confusion_matrix(y, y_pred)
report = classification_report(y, y_pred)


# In[6]:


print('x:', x, sep='\n')
print('y:', y, sep='\n', end='\n\n')
print('intercept:', model.intercept_)
print('coef:', model.coef_, end='\n\n')
print('p_pred:', p_pred, sep='\n', end='\n\n')
print('y_pred:', y_pred, end='\n\n')
print('score_:', score_, end='\n\n')
print('conf_m:', conf_m, sep='\n', end='\n\n')
print('report:', report, sep='\n')


# In[ ]:




