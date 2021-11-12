"""Helper functions for mapper analysis of plant RNAseq data.

This file contains helper functions used in the mapper analysis of the
plant RNAseq data.

Functions
---------
    - loaddata() : Loads and preprocesses factors and RNAseq data.
    - colorscale_from_matplotlib_cmap(): Generates a keppler mapper colorscale
        from matplotlib cmap.

Required packages
-----------------
    - numpy
    - pandas
    - scikit-learn

"""
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale


def loaddata(factor_file: str, rna_file: str, factor_list: list,
             logscale: bool = False, zscore: bool = True, group: bool = False) -> (object, list):
    """Load and Preprocess Factor and RNAseq Data.

    This function is used to load preprocess the data provided in two input
    files. First file contains factors data. Second file contains RNAseq data.
    The data points in two files are matched by a unique 'sra'.

    The function performs an inner join of two files on the 'sra'.
    This extracts records for which factors and RNAseq are both available.

    Then, preprocessing steps (scaling, normalization, grouping) are performed.
    A single dataframe, containing both factors and RNAseq is returned along
    with a list of orthogroup names.

    Parameters
    ----------
    factor_file : str
        full path to the csv file containing factor data
    rna_file : str
        full path to the csv file containint RNAseq data
    factor_list : list
        List of factor names to be kept (column names in factor_file)
    logscale : bool, optional
        Scale RNAseq data on a log2 scale (y = np.log2(x+1)).
        The default is True.
    zscore : bool, optional
        Standardize RNAseq data.
        Each 'sra' is normalized to 0 mean 1 standard deviation.
        The default is True.
    group : bool, optional
        Group records by sample_id.
        There's a one-to-many mapping between sample_id and sra.
        If true, factors data is grouped by sample id (keeping first record)
        RNAseq data is averaged across sra.
        The default is True.

    Returns
    -------
    df : pandas dataframe
        Contains factors and orthogroups data after preprocessing.
    ortho_names : list
        Contains list of names of orthogroups.
    """
    # Read the two files into pandas dataframes
    try:
        factor_input = pd.read_csv(factor_file)
    except FileNotFoundError:
        print('Could not find {}'.format(factor_file))
        return None, None
    else:
        fnr, fnc = factor_input.shape
        print('Factor input data shape: ({}, {})'.format(fnr, fnc))

    try:
        rna_input = pd.read_csv(rna_file)
    except FileNotFoundError:
        print('Could not find {}'.format(rna_file))
        return None, None
    else:
        rnr, rnc = rna_input.shape
        print('RNA input data shape: ({}, {})'.format(rnr, rnc))

    # Transpose the RNAseq data so that rows are indexed by 'sra'
    rna_input = rna_input.set_index('Orthogroup').T
    rna_input = rna_input.rename_axis('sra')\
        .rename_axis(None, axis='columns').reset_index()

    # Factor names and Orthogroup names
    ortho_names = (rna_input.columns[1:]).tolist()

    # Merge (inner join) factors and RNAseq data by "sra".
    df = pd.merge(factor_input, rna_input, on='sra', how='inner')

    # Print out the after-merge shapes of dataframes
    print('Dataframe shape after merge: {}'.format(df.shape))

    # Transform and standardize RNAseq data (z-score by row)
    if logscale:
        df.iloc[:, fnc:] = (df.iloc[:, fnc:]).apply(lambda x: np.log2(x+1.0))

    if zscore:
        df.iloc[:, fnc:] = scale(df.iloc[:, fnc:], axis=1)

    if group:
        df_grouped = df.groupby('sample_id')
        df = pd.concat([df_grouped[factor_list].first(),
                        df_grouped[ortho_names].mean()], axis=1)

        # Print out the shapes of dataframe after grouping
        print('Data shape after grouping: {}'.format(df.shape))

    return df, ortho_names


def colorscale_from_matplotlib_cmap(cmap):
    """Generate colorscale for keppler mapper from a matplotlib cmap.

    Parameters
    ----------
    cmap : matplotlib cmap
        For example matplotlib.pyplot.cm.plasma

    Returns
    -------
    colorscale : list of tuples
        range of values from 0 to 1 inclusive zipped with RGB strings.
        colorscale format compatible with keppler mapper.

    """
    cmap_list = [cmap(el) for el in np.arange(cmap.N)]
    rgb_strings = ["rgb({}, {}, {})".format(
            int(255 * el[0]), int(255 * el[1]), int(255 * el[2])
        ) for el in cmap_list
    ]

    idx_scale = np.linspace(0, 1., num=cmap.N)

    return list(zip(idx_scale, rgb_strings))


if __name__ == "__main__":

    datadir = '../data/SVA_Corrected_20210928'
    factorfile = datadir + '/clean_metadata.csv'
    rnafile = datadir + '/clean_RNAseq_sv_corrected.csv'
    factors = ['stress', 'tissue', 'family']
    levels = ['healthy', 'leaf', 'Poaceae']

    data, orthos = loaddata(factorfile, rnafile, factors)
    print(data.shape)
