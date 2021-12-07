#Adding all libraries to create a RF model
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#Import CSV for RF model
data = pd.read_csv('MapRateFiltered_tiny_75_seedling_deleted.csv')

#Build model
#distinguish between labels and features
features = data.iloc[:, 12:1001].values
labels = data.iloc[:, 8].values

#split the data
from sklearn.model_selection import train_test_split
f_train, f_test, l_train, l_test = train_test_split(features, labels, test_size = 0.25, random_state = 0)

#build the RF model
from sklearn.ensemble import RandomForestClassifier
classifier = RandomForestClassifier(n_estimators = 20, criterion = 'entropy', random_state = 0)
classifier.fit(f_train,l_train)
l_pred = classifier.predict(f_test)

from sklearn import metrics
print("Accuracy:",metrics.accuracy_score(l_test, l_pred))
