"""Factor Specific Genomic Analysis (FSGA).

In this script, we implement factor specific genomic analysis (FSGA) to plant
data. This follows the method from Nicolau et al.
Input data is in two files: First file contains several different factors
per sra / sample_id. Second file contains the RNAseq data.


Required packages
-----------------
    - numpy
    - scikit-learn

External scripts
-----------------
    - helper_functions.py

"""
import numpy as np
from sklearn.decomposition import PCA


def fsga_transform(data, ortho_names, factor, level):
    """FSGA transform RNA data w.r.t. specified factor and level.

    This function splits the samples into training and testing sets.
    The training set consists of samples that have the specified factor level.
    It then builds a linear model for the training set using PCA.
    First 'k' principal components that explain 90% of the variance in the
    training data are used.
    Then, the model is used to predict the testing set.
    Lastly, the residual vectors for training and testing sets are returned.

    Parameters
    ----------
    data : pandas dataframe
        Contains the factor and RNAseq data.
    ortho_names : list
        List of orthogroup names (column names in the dataframe)
    factor : str
        The factor w.r.t which to build the linear model.
    level : str
        The factor level used to determine the training set.

    Returns
    -------
    rna_residuals : pandas dataframe
        Contains residual vectors for both training and testing set
    idx_tr : list
        Contains indices for the training set
    idx_te : list
        Contains indices for the test set

    """
    rna_residuals = data[ortho_names].copy(deep=True)
    rna_residuals.reset_index(drop=True, inplace=True)
    idx_tr = (rna_residuals.index[data[factor] == level]).tolist()
    idx_te = (rna_residuals.index[data[factor] != level]).tolist()

    X = rna_residuals.iloc[idx_tr, :].to_numpy()
    Y = rna_residuals.iloc[idx_te, :].to_numpy()

    pca = PCA(n_components=0.9).fit(X)
    P = pca.components_.T @ pca.components_

    Xproj = X @ P
    X_residuals = X - Xproj

    Yproj = Y @ P
    Y_residuals = Y - Yproj

    rna_residuals.iloc[idx_tr, :] = X_residuals
    rna_residuals.iloc[idx_te, :] = Y_residuals

    return rna_residuals, idx_tr, idx_te


if __name__ == "__main__":

    datadir = '../data/SVA_Corrected_20210928'
    factorfile = datadir + '/clean_metadata.csv'
    rnafile = datadir + '/clean_RNAseq_sv_corrected.csv'
    factors = ['stress', 'tissue', 'family']
    levels = ['healthy', 'leaf', 'Poaceae']
    df, orthos = loaddata(factorfile, rnafile, factors)

    residuals, tr_idx, te_idx = fsga_transform(df, orthos, 'stress', 'healthy')

    print(np.linalg.norm(df[orthos], 'fro'))
    print(np.linalg.norm(residuals, 'fro'))
    print(np.linalg.norm(residuals.iloc[tr_idx, :], 'fro'))
    print(np.linalg.norm(residuals.iloc[te_idx, :], 'fro'))
