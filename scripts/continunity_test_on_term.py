import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from itertools import product

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

import warnings
warnings.filterwarnings('ignore')

if 'snakemake' in globals():
    from types import SimpleNamespace
    args = dict(input=snakemake.input[0],
                threshold=snakemake.params[0])
    args = SimpleNamespace(**args)
else:
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Get Stuff')

    parser.add_argument('--input',
                        type=str,
                        required=True,
                        help='Input Filename from Metaneighbor')
    parser.add_argument(
        '--threshold',
        type=float,
        default=.7,
        help='Threshold for binarizing the metaneighbor results (default : .7')
    args = parser.parse_args()


def isotonic_regression(dist):
    reg = IsotonicRegression(increasing='False').fit(dist['c'].values,
                                                     dist['r'].values)
    dist['fitC'] = reg.transform(dist['c'].values)
    dist['resid'] = np.abs(dist['r'].values - dist['fitC'].values)

    dist['max_resid'] = np.max(np.row_stack(
        [dist['fitC'], n_bins - dist['fitC']]),
                               axis=0)

    dist['norm_resid'] = dist['resid'].values / dist['max_resid'].values
    return dist

def get_binary_coord(df):
    return pd.DataFrame(np.where(df), index=['r', 'c']).T

def test_isotonicRegression(mat,
                            dist_func=get_binary_coord,
                            n_bins=16,
                            ret_val='norm_resid'):
    dist = isotonic_regression(dist_func(mat))
    return dist[ret_val].mean()


def is_monotonic_col(mat):
    return np.apply_along_axis(lambda x: len([i[0] for i in groupby(x)]), 0,
                               mat) <= 3


def percent_monotonic_cols(mat):
    monotonic = is_monotonic_col(mat)
    return monotonic.sum() / monotonic.shape[0]

def read_mn_res(fn):
    df = pd.read_csv(fn,
                     index_col=['Lineage', 'Study',
                                'Pseudotime Bin']).drop(columns=['Unnamed: 0'])
    df.columns = df.index
    return df


def prep_results(df):
    datasets = df.columns.levels[1]
    lineages = df.columns.levels[0]

    logging.info(
        f'Calculating results from # of Datasets : {datasets.shape[0]} ')
    logging.info(
        f'Calculating results from # of Lineages : {lineages.shape[0]} ')

    dataset_pairs = list(
        filter(lambda t: t[0] < t[1], list(product(datasets, datasets))))

    groups = list(
        filter(lambda t: t[1] < t[2],
               list(product(lineages, datasets, datasets))))

    results = pd.DataFrame(index=pd.MultiIndex.from_tuples(
        groups, names=['Lineage', 'Dataset1', 'Dataset2']),
                           columns=['Continuous', 'Monotonic'],
                           dtype=float)
    return results


def compute_pairwise_scores(mn_res, threshold=.7):
    results = prep_results(mn_res)
    for row in results.itertuples():
        lin, ds1, ds2 = row.Index
        mat = mn_res.loc[(lin, ds1), (lin, ds2)].values > threshold
        try:
            res_c = test_isotonicRegression(mat, dist_func=get_binary_coord)
            res_m = percent_monotonic_cols(mat)
        except:
            res_c = np.nan
            res_m = np.nan
        results.loc[row.Index, 'Continuous'] = res_c
        results.loc[row.Index, 'Monotonic'] = res_m
    return results


metaneighbor_output_fn = args.input
threshold = args.threshold


logging.info(f'Reading results from: {metaneighbor_output_fn.split("/")[-1]}')
mn_output = read_mn_res(metaneighbor_output_fn)
res = compute_pairwise_scores(mn_output, threshold=threshold)

output_fn = metaneighbor_output_fn.replace('MNUS', 'CONTINUITY')
logging.info(f'Saving Results to : {metaneighbor_output_fn.split("/")[-1]}')
res.to_csv(output_fn)
