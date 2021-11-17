## Python script to use Naive Bayes on the 3 gene test data

## import required packages
import pandas as pd
import numpy as np
import random

# Import Gaussian Naive Bayes classifier:
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score

## read in the data setting the first column as the index
# for 3 genes and 100 genes
dat = pd.read_csv("HRT841_3genes_test.csv", index_col=0)
dat100 = pd.read_csv("subset_tissue.csv", index_col=0)

# subset just the genes for the x data
xdata = dat[["AT1G22630","AT1G22620","AT1G22610"]]
xdata2 = dat100.iloc[:, 1:101]

# check the distributions of the expression data (Gaussian assumes Normal distribution)
pd.DataFrame.hist(data = xdata2, grid = False, bins = 40, figsize = (20,20))
# so definitely not all normally distributed but that is alright for now
# Categorical NB could be a good choice if the features are categorical


# now get the y data that I want
# will want to do more work here and might need to clean up a bit before my model is informative
y_tissues = dat["Tissue"]
y_stress = dat["Treatment"]
# really just care wild type or mutant here?
y_genotype = dat['Genotype']

y_tissues = dat100["Tissue"]
y_stress = dat100["Treatment"]
y_genotype = dat100['Genotype']

# features is the x data. labels is the y data.
# make test and train
# starting with tissue as the label
train, test, train_tissue, test_tissue = train_test_split(xdata2, y_tissues, test_size=0.33, random_state=42)

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

#potential issue that I am using Gaussian but I doubt we have normally distributed data
