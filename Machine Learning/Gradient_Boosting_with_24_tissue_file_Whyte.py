#!/usr/bin/env python
# coding: utf-8

# Gradient boosting is an ensamble method of supervised machine learning (ML) that consolidates the predictions of multiple trained models. 
# 
# Boosting refers to any model that bolsters low accuracy predictors by combining them to form stronger, more accurate predictors. Boosting models sequentially train its models so that each iteration builds upon the its predecessor, attempting to correct the errors that occured in the previous model. 
# 
# Gradient boosting functions by training models to correct the previous model's residual errors. The model starts with calculating a 'dummy estimator' by averaging the mean target value and makes predictions. The first residual is calculated by finding the difference between the predicted and actual value. 
# The next model attempts to predict the residuals of the first predictor (generally a decision tree with specific constraints). For each iteration, the base estimator's value is added to the predicted value of the decision tree to make a new prediction. The residuals are calculated again from the new prediction and the actual value. This process repeats until a specified threshold is met, or the residuals become very small.
# 
# This ML method can be used for classification or regression. Gradient boosting is useful for decreasing bias error and often outperforms Random Forest. However, it can overfit if your data has a lot of noise. Due to how the model functions, this method requires more processing power than other ML methods.

# In[ ]:


#Necessary imports
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.datasets import make_classification
from numpy import mean
from numpy import std
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import KFold
#(Optional) Necessary import if you encounter an error: "convert strings to floats"
from sklearn.preprocessing import LabelEncoder
urlfile="https://raw.githubusercontent.com/DanChitwood/PlantsAndPython_2021/main/Machine%20Learning/result_table_100_24_tissues.csv"

df = pd.read_csv(urlfile)
#Strings to float error fixing: Extract Tissue column, then re-label the column as strings
labels = df["Tissue"].values
label_encoder = LabelEncoder().fit(np.unique(labels))
df['Tissue'] = label_encoder.transform(labels)


# In[ ]:


#(Optional) Check the column headings for indexing to change the number of genes you would 
# like to use or the category you're predicting
df.head()


# In[ ]:


#assign X and Y variables. X must include what you are predicting, Y will just be the values to predict
X = df.iloc[:, 2:102].values
y = df.iloc[:, 102].values

#Define features
kf = KFold(n_splits=5,random_state=42,shuffle=True)
train_percent = 0.75
# create training and testing subsets
train = np.random.choice(a=[True, False], size=len(y), p=[train_percent, 1-train_percent])
test = np.invert(train)
for train_index,val_index in kf.split(X):
    X_train,X_test = X[train],X[test],
    y_train,y_test = y[train],y[test]


# In[ ]:


#(Optional) this statement retrieves the adjustable parameters in addition to revealing the default settings
#The following parameters may be the most useful to change:
#     Criterion: The loss function used to find the optimal feature and threshold to split the data
#     learning_rate: this parameter scales the contribution of each tree
#     max_depth: the maximum depth of each tree
#     n_estimators: the number of trees to construct
#     init: the initial estimator. By default, it is the log(odds) converted to a probability
gradient_booster.get_params()

#Gradient booster model
gradient_booster = GradientBoostingClassifier(learning_rate=0.1, max_depth=20, random_state=10)


# In[ ]:


#Fit the model
GB_fit=gradient_booster.fit(X_train,y_train)
#Print the accuracy of the model
print(classification_report(y_test,gradient_booster.predict(X_test)))


# In[ ]:


#Visualizing the model
plt.figure(figsize=(12,8)) 
tree.plot_tree(GB_fit, filled=True, fontsize=10)
plt.show()

