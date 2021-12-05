from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler 
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import datasets # import dataset where is iris
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import load_iris
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, normalize


dataset = pd.read_csv('result_table_100_24_tissues.csv')  


# Use head() function to return the first 5 rows: 
dataset.head()


X = dataset.iloc[:, 2:102].values
y = dataset.iloc[:, 102].values

X.shape


X_norm = normalize(X)
X = pd.DataFrame(X_norm)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)

scaler = StandardScaler()
scaler.fit(X_train)


X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test) 


# Neural network
clf = MLPClassifier(solver='lbfgs', #optimizador lbfgs
                    hidden_layer_sizes=(10,2), random_state=1)
clf.fit(X_train, y_train)
pred = clf.predict(X_test)
cnf_matrix = confusion_matrix(y_test, pred)
print('10 neurons 2 layer: ' + str(accuracy_score(y_test, pred)))
