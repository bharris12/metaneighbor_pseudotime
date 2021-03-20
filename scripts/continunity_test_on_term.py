"""Summary

Attributes
----------
beta : TYPE
    Description
binary_threshold : TYPE
    Description
metaneighbor_output_fn : TYPE
    Description
mn_output : TYPE
    Description
output_fn : TYPE
    Description
res : TYPE
    Description
soft_threshold : TYPE
    Description
metaneighbor_output_fn : TYPE
Description
mn_output : TYPE
Description
output_fn : TYPE
Description
res : TYPE
Description
threshold : TYPE
Description

Deleted Attributes
------------------
threshold : TYPE
    Description
"""
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from itertools import product, groupby

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

import warnings
warnings.filterwarnings('ignore')

if 'snakemake' in globals():
    from types import SimpleNamespace
    args = dict(input=snakemake.input[0],
                binary_threshold=snakemake.params[0],
                outdir=snakemake.params[1],
                dist_func=snakemake.params[2],
                beta=snakemake.params[3],
                soft_threshold=snakemake.params[4])
    args = SimpleNamespace(**args)
else:
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Get Stuff')

    parser.add_argument('--input',
                        type=str,
                        required=True,
                        help='Input Filename from Metaneighbor')
    parser.add_argument(
        '--binary-threshold',
        type=float,
        default=.7,
        help='Threshold for binarizing the metaneighbor results (default : .7')
    parser.add_argument('--outdir', type=str, help='Output directory')
    parser.add_argument('--dist-func',
                        type=str,
                        default='weighted_dist',
                        help='Function for using to ')
    parser.add_argument('--beta',
                        type=int,
                        default=1,
                        help='Beta for raising weighted dist AUROC to')
    parser.add_argument(
        '--soft-threshold',
        type=float,
        default=.5,
        help=
        'Threshold for selecting coordinates to keep for scaling before fitting weithed regression'
    )
    args = parser.parse_args()


def read_mn_res(fn):
    """Summary
    
    Parameters
    ----------
    fn : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    df = pd.read_csv(fn,
                     index_col=['Lineage', 'Study',
                                'Pseudotime Bin']).drop(columns=['Unnamed: 0'])
    df.columns = df.index
    return df


def binary_dist(df):
    """Get coordinates for True values in df
    
    
    Parameters
    ----------
    df : TYPE
        Description
    df {pd.Dataframe} -- Pandas DataFrame with boolean values
    
    Returns
    -------
    pd.DataFrame -- Dataframe with r and c columns corresponding to the locations of the true values
    """
    return pd.DataFrame(np.where(df), index=['r', 'c']).T


def weighted_dist(df, beta, soft_threshold):
    if type(df) is not pd.DataFrame:
        df = pd.DataFrame(df)

    df.index.name = 'r'
    df.columns.name = 'c'
    df = df.reset_index().melt(id_vars='r', var_name='c', value_name='weight')
    df = df[df.weight > soft_threshold]
    df.weight = df.weight**beta
    return df


def isotonic_regression(dist, n_bins, x='c', y='r'):
    reg = IsotonicRegression(increasing='False').fit(
        dist[x].values,
        dist[y].values,
        sample_weight=dist['weight'].values
        if 'weight' in dist.columns else None)
    dist.loc[:, 'fitC'] = reg.transform(dist[x].values)
    dist.loc[:, 'resid'] = np.abs(dist[y].values - dist['fitC'].values)

    max_test = np.vectorize(lambda x, y: np.max([x, y]))
    dist.loc[:, 'max_resid'] = max_test(dist['fitC'].values,
                                        n_bins - dist['fitC'].values)
    dist.loc[:, 'norm_resid'] = dist['resid'].values / dist['max_resid'].values
    return dist


def ITregression(mat, n_bins, dist_func, beta, binary_threshold,
                 soft_threshold):
    if dist_func == 'weighted_dist':
        dist = weighted_dist(mat, beta, soft_threshold)
    elif dist_func == 'binary_dist':
        dist = binary_dist(mat > binary_threshold)
    else:
        raise NotImplementedError(dist_func)
    dist = isotonic_regression(dist, n_bins)
    return dist['norm_resid'].mean()


def is_monotonic_col(mat):
    """Summary
    
    Parameters
    ----------
    mat : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    return np.apply_along_axis(lambda x: len([i[0] for i in groupby(x)]), 0,
                               mat) <= 3


def percent_monotonic_cols(mat):
    """Summary
    
    Parameters
    ----------
    mat : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    monotonic = is_monotonic_col(mat)
    return monotonic.sum() / monotonic.shape[0]


def prep_results(df):
    """Summary
    
    Parameters
    ----------
    df : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    datasets = df.columns.levels[1]
    lineages = df.columns.levels[0]

    logging.info(
        f'Calculating results from # of Datasets : {datasets.shape[0]} ')
    logging.info(
        f'Calculating results from # of Lineages : {lineages.shape[0]} ')

    groups = list(
        filter(lambda t: t[1] < t[2],
               list(product(lineages, datasets, datasets))))

    results = pd.DataFrame(index=pd.MultiIndex.from_tuples(
        groups, names=['Lineage', 'Dataset1', 'Dataset2']),
                           columns=['Continuous', 'Monotonic'],
                           dtype=float)
    return results


def compute_pairwise_scores(mn_res, dist_func, beta, binary_threshold,
                            soft_threshold):
    """Summary
    
    Parameters
    ----------
    mn_res : TYPE
        Description
    dist_func : TYPE
        Description
    beta : TYPE
        Description
    binary_threshold : TYPE
        Description
    soft_threshold : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    
    Deleted Parameters
    ------------------
    threshold : float, optional
        Description
    """
    results = prep_results(mn_res)
    for row in results.itertuples():
        lin, ds1, ds2 = row.Index
        mat = mn_res.loc[(lin, ds1), (lin, ds2)].values
        n_bins = mat.shape[0]
        try:
            res_c = ITregression(mat, n_bins, dist_func, beta,
                                 binary_threshold, soft_threshold)
        except:
            res_c = np.nan
        try:
            res_m = percent_monotonic_cols(mat > binary_threshold)
        except:
            res_m = np.nan

        results.loc[row.Index, 'Continuous'] = res_c
        results.loc[row.Index, 'Monotonic'] = res_m
    return results


metaneighbor_output_fn = args.input
binary_threshold = args.binary_threshold
soft_threshold = args.soft_threshold
beta = args.beta
dist_func = args.dist_func

logging.info(f'Reading results from: {metaneighbor_output_fn.split("/")[-1]}')
mn_output = read_mn_res(metaneighbor_output_fn)

res = compute_pairwise_scores(mn_output, dist_func, beta, binary_threshold,
                              soft_threshold)
output_fn = args.outdir + metaneighbor_output_fn.replace(
    'MNUS', 'CONTINUITY').split('/')[-1]
logging.info(f'Saving Results to : {output_fn}')
res.to_csv(output_fn)
