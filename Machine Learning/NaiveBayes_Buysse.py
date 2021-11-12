## Python script to use Naive Bayes on the 3 gene test data

## import required packages
import pandas as pd
import numpy as np
import random

# Import Gaussian Naive Bayes classifier:
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn.naive_bayes import CategoricalNB

## read in the data setting the first column as the index
dat = pd.read_csv("HRT841_3genes_test.csv", index_col=0)
print(dat.head())
# subset just the genes for the x data
xdata = dat[["AT1G22630","AT1G22620","AT1G22610"]]
xdata.head()

# now get the y data that I want
# will want to do more work here and might need to clean up a bit before my model is informative
y_tissues = dat["Tissue"]
y_stress = dat["Treatment"]
# really just care wild type or mutant here?
y_genotype = dat['Genotype']

# features is the x data. labels is the y data.
# make test and train
# starting with tissue as the label
train, test, train_tissue, test_tissue = train_test_split(xdata, y_tissues, test_size=0.33, random_state=42)

# Initialize classifier:
gnb = GaussianNB()
cnb = CategoricalNB()

# Train the classifier:
model = gnb.fit(train, train_tissue)
model2 = cnb.fit(train, train_tissue)
# Make predictions with the classifier:
predictive_tissues = gnb.predict(test)
pred_2 = cnb.predict(test)
print(predictive_tissues)

# Evaluate label (subsets) accuracy:
print(accuracy_score(test_tissue, predictive_tissues))

print("GaussianNB: Number of mislabeled points out of a total %d points : %d"
      % (test.shape[0], (test_tissue != predictive_tissues).sum()))

print(accuracy_score(test_tissue, pred_2))

print("CateboricalNB: Number of mislabeled points out of a total %d points : %d"
      % (test.shape[0], (test_tissue != pred_2).sum()))

#potential issue that I am using Gaussian but I doubt we have normally distributed data
