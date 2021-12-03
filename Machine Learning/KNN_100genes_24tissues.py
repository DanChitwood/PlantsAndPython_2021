#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier, NeighborhoodComponentsAnalysis
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix


# In[3]:


dataset24tissues = pd.read_excel('result_table_100_24_tissues.xlsx') # data set
 

X = dataset24tissues.iloc[:, 2:102].values
y = dataset24tissues.iloc[:, 102].values


# In[5]:


# Assign values to the X and y variables:
X = dataset24tissues.iloc[:, 2:102].values
y = dataset24tissues.iloc[:, 102].values


# In[6]:


# Split dataset into random train and test subsets:
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)

#train_test_split is a function in Sklearn model selection for splitting data arrays into two subsets


# In[7]:


n_neighbors = 5 # K = 5

# Reduce dimension to 5 with PCA
pca = make_pipeline(StandardScaler(), PCA(n_components=5))

# Reduce dimension to 5 with LinearDiscriminantAnalysis
lda = make_pipeline(StandardScaler(), LinearDiscriminantAnalysis(n_components=5))

# Reduce dimension to 5 with NeighborhoodComponentAnalysis
nca = make_pipeline(
    StandardScaler(),
    NeighborhoodComponentsAnalysis(n_components=5),
)

# Use a nearest neighbor classifier to evaluate the methods
knn = KNeighborsClassifier(n_neighbors=n_neighbors)

# Make a list of the methods to be compared
dim_reduction_methods = [("PCA", pca), ("LDA", lda), ("NCA", nca)]


# In[8]:


knn.fit(X_train, y_train) # fit the model


# In[9]:


y_predict = knn.predict(X_test) # predict the testing population


# In[10]:


print(classification_report(y_test, y_predict)) #print the prescitions


# In[11]:


# this graphic is not vey useful, I think that is because is in 2 dimentions

for i, (name, model) in enumerate(dim_reduction_methods):
    plt.figure()

    # Fit the method's model
    model.fit(X_train, y_train)

    # Fit a nearest neighbor classifier on the embedded training set
    #knn.fit(model.transform(X_train), y_train)
    knn.fit(X_train, y_train)
    # Compute the nearest neighbor accuracy on the embedded test set
    #acc_knn = knn.score(model.transform(X_test), y_test)
    acc_knn = knn.score(X_test, y_test)  


    # Embed the data set in 2 dimensions using the fitted model
    #X_embedded = model.transform(X)
    X_embedded = (X)

    # Plot the projected points and show the evaluation score
    sns.relplot(X_embedded[:, 0], X_embedded[:, 1],
              hue=y, kind="scatter", palette=sns.hls_palette(24, l=0.7, s=1))
    plt.title(
        "{}, KNN (k={})\nTest accuracy = {:.2f}".format(name, n_neighbors, acc_knn)
    )
plt.show()

