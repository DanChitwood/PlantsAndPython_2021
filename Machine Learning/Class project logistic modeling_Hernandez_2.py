#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
import pandas as pd


# In[2]:


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


# In[3]:


# model = LogisticRegression(solver='liblinear', random_state=0).fit(x, y)
model = LogisticRegression(solver='liblinear', random_state=0)


# In[4]:


model.fit(x, y)


# In[5]:


model.classes_


# In[6]:


model.intercept_


# In[7]:


model.coef_


# In[8]:


#evaluate the model
model.predict_proba(x)


# In[9]:


model.predict(x)


# In[10]:


model.score(x, y)


# In[11]:


confusion_matrix(y, model.predict(x))


# In[12]:


cm = confusion_matrix(y, model.predict(x))

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(cm)
ax.grid(False)
ax.xaxis.set(ticks=(0, 1), ticklabels=('Predicted 0s', 'Predicted 1s'))
ax.yaxis.set(ticks=(0, 1), ticklabels=('Actual 0s', 'Actual 1s'))
ax.set_ylim(1.5, -0.5)
for i in range(2):
    for j in range(2):
        ax.text(j, i, cm[i, j], ha='center', va='center', color='red')
plt.show()


# In[13]:


print(classification_report(y, model.predict(x)))


# In[14]:


#improve the model
model = LogisticRegression(solver='liblinear', C=10.0, random_state=0)
model.fit(x, y)


# In[15]:


model.intercept_
model.coef_
model.predict_proba(x)
model.predict(x)


# In[17]:


model.score(x,y)
confusion_matrix(y, model.predict(x))
print(classification_report(y, model.predict(x)))


# In[ ]:




