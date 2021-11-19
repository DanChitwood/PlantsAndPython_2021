## Python script to use Naive Bayes on the 3 gene test data

## import required packages
import pandas as pd
import numpy as np
import random

# Import Gaussian Naive Bayes classifier:
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

## read in the data setting the first column as the index
# the 100 genes and 24 tissues dataset 9on github)
dat24 = pd.read_csv("result_table_100_24_tissues.csv", index_col=0)

# subset just the genes for the x data
xdata = dat24.iloc[:, 1:101]

# check the distributions of the expression data
# this takes up a lot of space so only do once to get an idea
#pd.DataFrame.hist(data = xdata, grid = False, bins = 40, figsize = (20,20))
# so definitely not all normally distributed but that is alright for now
# Categorical NB could be a good choice if the features are categorical (not what we have)


# now get the y data that I want
# will want to do more work here and might need to clean up a bit before my model is informative
unique_tissues = dat24["Tissue"].unique()
print(len(unique_tissues))
y_tissues = dat24["Tissue"]

# features is the x data. labels is the y data.
# make test and train
# starting with tissue as the label
train, test, train_tissue, test_tissue = train_test_split(xdata, y_tissues, test_size=0.25, random_state=42)

# Initialize classifier:
gnb = GaussianNB()

# Train the classifier:
model = gnb.fit(train, train_tissue)

# Make predictions with the classifier:
predictive_tissues = gnb.predict(test)
print(predictive_tissues)

# Evaluate label (subsets) accuracy:
print(accuracy_score(test_tissue, predictive_tissues))

print("GaussianNB: Number of mislabeled points out of a total %d points : %d"
      % (test.shape[0], (test_tissue != predictive_tissues).sum()))

# potential issue in choosing Gaussian. Do we have a prior expectation that the relationship of expression to tissue type would be normal?

##### from Paulo's code for testing performance #####
# In the field of machine learning and specifically the problem of statistical classification, 
# a confusion matrix, also known as an error matrix

# Each row of the matrix represents the instances in an actual class while each column 
# represents the instances in a predicted class,

print(confusion_matrix(test_tissue, predictive_tissues))


print(classification_report(test_tissue, predictive_tissues)) 
