## Python script to use Naive Bayes
## HAS NOT BEEN UPDATED TO USE THE TISSUE SUBSET FROM GROUP TWO
## this document uses a subset of 100 genes and 24 tissue types made by Harry

## Author: Sophie Buysse with some code from Paulo to look at accuracy

## import required packages
import pandas as pd
import numpy as np
import random

# Import Gaussian Naive Bayes classifier:
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

## read in the data and set the first column as the index
# this code uses the 100 genes and 24 tissues dataset (on github in the ML folder)
dat24 = pd.read_csv("result_table_100_24_tissues.csv", index_col=0)

# subset just the genes for the x data
# the genes are our features (or predictors) so they are the x data we will be using
xdata = dat24.iloc[:, 1:101]

# check the distributions of the expression data
# this takes up a lot of space so only do once to get an idea
# Gaussian NB assumes the features are normally distributed, so make a histogram for each feature to see if that is the case
#pd.DataFrame.hist(data = xdata, grid = False, bins = 40, figsize = (20,20))
# so definitely not all normally distributed but that is alright for now
# because we would be using a different set of features after feature selection

# Categorical NB could be a good choice if the features are categorical (not what we have)


## now get the y data that we want
# the y data is the classes we want to sort in to
# will want to do more work here and might need to clean up a bit before my model is informative
# depends on what a more final dataset looks like

# double check there are the correct number of tissues as expected
unique_tissues = dat24["Tissue"].unique()
print(len(unique_tissues))

#subset the tissue column to a new object
y_tissues = dat24["Tissue"]

# features is the x data. classes is the y data.
# make test and train
# starting with tissue as the label
# test should be smaller than the train size
train, test, train_tissue, test_tissue = train_test_split(xdata, y_tissues, test_size=0.25, random_state=42)

# Initialize classifier:
gnb = GaussianNB()

# Train the classifier:
model = gnb.fit(train, train_tissue)

# Make predictions with the classifier:
# (predict the tissue based on the gene expression)
# do this on the test dataset
predictive_tissues = gnb.predict(test)
print(predictive_tissues)

# Evaluate label (subsets) accuracy:
# (did it classify as expected?)
print(accuracy_score(test_tissue, predictive_tissues))

print("GaussianNB: Number of mislabeled points out of a total %d points : %d"
      % (test.shape[0], (test_tissue != predictive_tissues).sum()))

##### from Paulo's code for testing performance #####
# In the field of machine learning and specifically the problem of statistical classification, 
# a confusion matrix, also known as an error matrix

# Each row of the matrix represents the instances in an actual class while each column 
# represents the instances in a predicted class,

print(confusion_matrix(test_tissue, predictive_tissues))


print(classification_report(test_tissue, predictive_tissues)) 


### final notes on how this went:
# not a great model fit; a different model is probably better
# could take more time to think about a prior and if Gaussian is even appropriate
