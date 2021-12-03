"""
Fit logistic regression using L1, L2, and Elastic Net penalties.

Penalty strength is determined via cross-validation
"""
import warnings
import numpy as np
import pandas as pd 
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder 
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.exceptions import ConvergenceWarning

# Get rid of annoying warnings
warnings.filterwarnings('ignore', category=ConvergenceWarning)


def get_data(label="Tissue"):
    """
    Return training and test data.

    Do all cleaning, splitting, and standardization here.

    Most here lifted from Paolo's file

    Parameters:
    -----------
        label: str
            What we are trying to predict. "Tissue" by default

    Returns:
    --------
    tuple
        X_train, Y_
    """
    # df = pd.read_excel('subset_tissue.xlsx', engine='openpyxl')
    df = pd.read_csv("result_table_100_24_tissues.csv")

    print("\nLabel Distribution:\n-------------------")
    print(df['Tissue'].value_counts(), "\n")

    print("Total number of Samples:", df.shape[0], "\n")

    features = df.iloc[:, 2:102].values
    labels = df[label].values

    # Some methods expect numeric labels...
    label_encoder = LabelEncoder().fit(np.unique(labels))
    labels = label_encoder.transform(labels)

    X_train, X_test, y_train, y_test = train_test_split(features, labels, stratify=labels, test_size=0.20)

    # Only scale based on training data so no leakage
    scaler = StandardScaler()
    scaler.fit(X_train)

    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    return X_train, y_train, X_test, y_test



def eval_model(model, X_test, y_test):
    """
    Evaluate specific model using predictions on test set
    """
    preds = model.predict(X_test)

    print(confusion_matrix(preds, y_test))
    print(classification_report(preds, y_test)) 
    print("\n\n")


def fit_logistic(X_train, y_train, X_test, y_test, penalty="l2"):
    """
    """
    reg_strength = [1e-3, 1e-2, 1e-1, 1, 1e2, 1e3]

    print(f"Fitting logistic regression with {penalty} penalty...")

    model = LogisticRegressionCV(Cs=reg_strength, cv=5, max_iter=250, random_state=42).fit(X_train, y_train)

    eval_model(model, X_test, y_test)



def main():
    X_train, y_train, X_test, y_test = get_data()

    # for p in ['l1', 'l2', 'elasticnet']:
    #     fit_logistic(X_train, y_train, X_test, y_test, penalty=p)

    fit_logistic(X_train, y_train, X_test, y_test, penalty='l1')
    

if __name__ == "__main__":
    main()
