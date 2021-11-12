#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Import libraries and classes required for this example:
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler 
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import datasets # import dataset where is iris
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import load_iris


# In[3]:


dataset = pd.read_excel('subset_tissue.xlsx')  


# In[4]:


dataset


# In[5]:


# Use head() function to return the first 5 rows: 
dataset.head()


# In[19]:


# Assign values to the X and y variables:
X = dataset.iloc[:, 2:102].values
y = dataset.iloc[:, 102].values


# In[20]:


X.shape


# In[39]:


dataset.shape # dimensions


# In[21]:


# Split dataset into random train and test subsets:
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)

#train_test_split is a function in Sklearn model selection for splitting data arrays into two subsets


# In[22]:


print(X_train.shape)
print(X_test.shape)


# In[23]:


# Standardize features by removing mean and scaling to unit variance:
scaler = StandardScaler()
scaler.fit(X_train)


# In[24]:


X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test) 


# In[25]:


# Use the KNN classifier to fit data:
# K-Nearest Neighbors (KNN) – a simple classification algorithm, 
# where K refers to the square root of the number of training records.
## Scikit-Learn

classifier = KNeighborsClassifier(n_neighbors=5)
# ‘n_neighbors‘ are the number of neighbors that will vote for the class of the target point; 
# default number is 5. An odd number is preferred to avoid any tie.

classifier.fit(X_train, y_train) 


# In[26]:


# Predict y data with classifier: 
y_predict = classifier.predict(X_test)


# In[28]:


plt.plot( y_test, y_predict,  linestyle='', marker='o', markersize = 10)
plt.xlabel('Tissue observed')
plt.ylabel('Tissue predicted')
plt.show()


# In[29]:


ax = sns.stripplot(y_test, y_predict, jitter= 0.2, s = 9, marker="D", alpha=.5 )

ax.set(xlabel='Tissue', ylabel='Predictions')


# In[30]:


# In the field of machine learning and specifically the problem of statistical classification, 
# a confusion matrix, also known as an error matrix

# Each row of the matrix represents the instances in an actual class while each column 
# represents the instances in a predicted class,

print(confusion_matrix(y_test, y_predict))


# In[31]:


print(classification_report(y_test, y_predict)) 


# In[36]:


dataset.shape

