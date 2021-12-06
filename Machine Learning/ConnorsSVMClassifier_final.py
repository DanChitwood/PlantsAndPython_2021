#!/usr/bin/env python
# coding: utf-8

# In[31]:


#Importing needed packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.pipeline import Pipeline
from sklearn.metrics import plot_confusion_matrix
from sklearn.preprocessing import StandardScaler
get_ipython().run_line_magic('matplotlib', 'inline')


# In[32]:


#importing data
dataset = pd.read_excel('tissue15genes.xlsx')  
dataset.head()


# In[33]:


#setting X and Y values; X is all gene expression data, y is tissue type
X=dataset.iloc[:,3:102]
y = dataset["Tissue"]


# In[34]:


#scaling the data
scaler = StandardScaler()
model = scaler.fit(X)
scaled_data = model.transform(X)


# In[35]:


#Split the data into training and testing
X_train, X_test, y_train, y_test = train_test_split(scaled_data,y,test_size=0.2)


# In[36]:


#train the SVM model
svcclassifier=SVC(kernel='linear')
svcclassifier.fit(X_train,y_train.values.ravel())


# In[37]:


#Predict y values
y_pred=svcclassifier.predict(X_test)


# In[38]:


#Evaluate the model
print(confusion_matrix(y_test,y_pred))


# In[39]:


print(classification_report(y_test,y_pred))


# In[ ]:




