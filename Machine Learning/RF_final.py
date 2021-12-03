"""This code uses Random Forest Classification to predict what tissue type a specific Arabidopsis 
sample is based on gene expression data. There are 24 classes of labels, the instances are the 
samples, and the features are the genes. Each instance has data about gene expression for each
gene."""

#Importing Libraries
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
import joblib
#import for checking our RF model
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
#print('Libraries Imported') - check

#Load in the dataset
dataset = pd.read_csv('result_table_100_24_tissues.csv', header=0)
#confirm data read in correctly
#print('Shape of the dataset: ' + str(dataset.shape))
#print(dataset.head())

#variables
#number of features, label column, n_estimators, number of labels
features = 102
label_column = 102
n = 10
n_labels = 24

#Turn Tissue types into numbers
numbers = pd.factorize(dataset['Tissue'])
dataset.Tissue = numbers[0]
definitions = numbers[1]
#print(dataset.Tissue.head()) - check
#print(definitions) - check

#Splitting data into independent (features - genes) and dependent (labels - tissue types) variables
genes = dataset.iloc[:,2:features].values
tissue_types = dataset.iloc[:,label_column].values
#print('The independent features set: ')
#print(genes[:5,:]) - check
#print(tissue_types[:8]) - check

#Use 75% for training and 25% of data for testing
#random state assigns random distribution of data
genes_train, genes_test, tissue_types_train, tissue_types_test = train_test_split(genes, tissue_types, test_size=0.25, random_state=21)

#Feature scaling - subtracts mean value of observed and divides by unit variance of observed
scaler = StandardScaler()
genes_train = scaler.fit_transform(genes_train)
genes_test = scaler.transform(genes_test)

#Train the model
#Fit RF Classification to Training Set
#n_estimators is how deep the model goes
classifier = RandomForestClassifier(n_estimators= n, criterion='entropy', random_state= 42)
classifier.fit(genes_train, tissue_types_train)

#Evaluate performance
#Predict labels (tissue_types) of test data (genes_test) - use predict to do so
tissue_types_predict = classifier.predict(genes_test)
#Reverse factorize
reverse_factor = dict(zip(range(n_labels),definitions))
tissue_types_test = np.vectorize(reverse_factor.get)(tissue_types_test)
tissue_types_predict = np.vectorize(reverse_factor.get)(tissue_types_predict)
#Make a cross tabulation to see the data
print(pd.crosstab(tissue_types_test, tissue_types_predict, rownames=['Actual Tissue Type'], colnames=['Predicted Tissue Types']))
#Check accuracy
accuracy = accuracy_score(tissue_types_test, tissue_types_predict)
print(accuracy)

#Get and reshape confusion matrix data - option: put this into a heatmap to visualize 
matrix = confusion_matrix(tissue_types_test, tissue_types_predict)
matrix = matrix.astype('float')
matrix.sum(axis= 1) [:, np.newaxis]
#print(matrix)

"""Written by Brianna Brown, performance on dataset is 0.873 with n_estimators = 10"""