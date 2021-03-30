import numpy as np
import pandas as pd
import scanpy as sc
import pymn
from pathlib import Path
import mkl
import logging
from itertools import product
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

##Functions

get_study = np.vectorize(lambda x: x.split('|')[0])
get_lineage = np.vectorize(lambda x: x.split('|')[-1].split('_')[0])
get_bin = np.vectorize(lambda x: int(x.split('|')[-1].split('_')[-1]))
joining = np.vectorize(lambda x, y: f'{x}_{y}')


def compute_bins(pseudotime, n_bins):
    return pd.cut(pseudotime, n_bins,
                  labels=np.arange(n_bins).astype(str)).astype(str)


def create_multi_index(df):
    s = get_study(df.index)
    l = get_lineage(df.index)
    b = get_bin(df.index)
    expected_index = [
        f'{a}|{b}_{c}'
        for a, b, c in list(product(np.unique(s), np.unique(l), np.unique(b)))
    ]
    df = df.reindex(expected_index).reindex(expected_index, axis=1)
    s = get_study(df.index)
    l = get_lineage(df.index)
    b = get_bin(df.index)

    df.index = pd.MultiIndex.from_arrays(
        [l, s, b], names=['Lineage', 'Study', 'Pseudotime Bin'])
    df.columns = pd.MultiIndex.from_arrays(
        [l, s, b], names=['Lineage', 'Study', 'Pseudotime Bin'])
    return df


def output_mn_bin_lineage_res(df, fn):
    mn_combined = df.copy()

    mn_combined = create_multi_index(mn_combined)

    mn_combined = mn_combined.sort_index().sort_index(axis=1).droplevel(
        ['Lineage', 'Study'], axis=1)

    mn_combined.columns = np.arange(mn_combined.shape[1])

    mn_combined.reset_index(inplace=True)
    mn_combined.to_csv(fn)


if 'snakemake' in globals():
    from types import SimpleNamespace
    args = dict(adata_fn=snakemake.input[0],
                term=snakemake.input[1],
                join_bin=snakemake.params[0],
                rank_by=snakemake.params[1],
                dataset_name=snakemake.params[2],
                n_bins=snakemake.params[3],
                output_path=snakemake.params[4],
                study_col=snakemake.params[5],
                pseudotime_col=snakemake.params[6],
                max_genes=snakemake.params[7],
                slow=snakemake.params[8],
                term_name=snakemake.wildcards[0],
                n_threads=snakemake.threads)
    args = SimpleNamespace(**args)
else:
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Get Stuff')

    parser.add_argument('--n-bins', type=int, help='Number of bins to use')
    parser.add_argument('--adata-fn',
                        type=str,
                        required=True,
                        help='Filename of Dataset')
    parser.add_argument('--dataset-name',
                        type=str,
                        required=True,
                        help='Name of Dataset')
    parser.add_argument(
        '--term',
        type=str,
        help=
        'Name of GO term or Gene list to compute MN on or filename with genes',
        default='HVG')
    parser.add_argument(
        '--term-name',
        type=str,
        help='If passing file with gene list, must provide name',
        default=None)
    parser.add_argument('--pseudotime-col',
                        type=str,
                        default='monocle_ps',
                        help='Column in adata.obs with Pseudotime data')
    parser.add_argument('--n-threads',
                        default=8,
                        type=int,
                        help='Number of threads to use')
    parser.add_argument('--study-col',
                        default='study_id',
                        help='Col in adata.obs that has the study name in it')
    parser.add_argument(
        '--output-path',
        default='/home/bharris/pseudotime/data/monocle/metaneighbor/')
    parser.add_argument(
        '--join-bin',
        type=str,
        help='Column to Join Bin to for additionally CT delineation')
    parser.add_argument(
        '--rank-by',
        type=str,
        nargs='+',
        help='Columns in dataframe to groupby for ranking pseudotime')
    parser.add_argument(
        '--slow',
        action='store_false',
        default=True,
        help=
        'Whether to run in Original Slow mode, default fast approximate version'
    )
    parser.add_argument('--max-genes',
                        default=100,
                        help='Max number of genes in the term')
    args = parser.parse_args()

adata_fn = args.adata_fn
dataset_name = args.dataset_name
n_bins = args.n_bins
study_col = args.study_col
n_threads = args.n_threads
study_col = args.study_col
output_path = args.output_path
pseudotime_col = args.pseudotime_col
join_bin = args.join_bin
rank_by = args.rank_by
term = args.term
term_name = args.term_name
max_genes = args.max_genes
slow = args.slow

mkl.set_num_threads(n_threads)
logging.info('Loading Dataset')
adata = sc.read(adata_fn)
adata = adata[adata.obs[join_bin] !='nan']
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
bin_col = f'ps_bin_{n_bins}'
if rank_by is not None:
    logging.info(f'Ranking pseudotime by {rank_by}')
    adata.obs[pseudotime_col] = adata.obs.groupby(
        rank_by)[pseudotime_col].rank(pct=True)

adata.obs[bin_col] = compute_bins(adata.obs[pseudotime_col], n_bins)

## Joining bins to a categorical column in the .obs data
if join_bin is not None:
    assert join_bin in adata.obs.columns, 'Join bin must be in the columns of the data'
    logging.info(f'Joining Bins to {join_bin}')
    adata.obs[bin_col] = joining(adata.obs[join_bin].values,
                                 adata.obs[bin_col].values)
else:
    adata.obs[bin_col] = joining([dataset_name] * adata.shape[0],
                                 adata.obs[bin_col].values)

## Get gene list genes
if Path(term).is_file():
    assert args.term_name != None, 'Must pass term name when passing file for term'
    term_genes = np.genfromtxt(term, dtype=str)
    term = args.term_name
elif term != 'HVG':
    logging.info('Reading in Gene Ontology')
    go_mouse = pd.read_hdf('/home/bharris/GO_data/go_mouse_nw.hdf5', 'go')
    assert term in go_mouse.columns, 'Term is not in GO, and HVG not provided'
    term_genes = go_mouse.index[go_mouse[term].astype(bool)].values
else:
    logging.info('Computing Highly Variable Genes')
    term_genes = adata.var_names[pymn.variableGenes(adata,
                                                    study_col,
                                                    return_vect=True)]

assert term_genes.shape[0] > 10, 'Must use a term with at least 10 genes'
assert term_genes.shape[
    0] < max_genes or term == 'HVG', f'Must use a term with less than {max_genes}'

adata.var['highly_variable'] = np.in1d(adata.var_names, term_genes)
logging.info(f'Running MNUS on {term} of size {term_genes.shape[0]}')
adata.obs[study_col] = adata.obs[study_col].astype(str)

##Running MetaNeighborUS on gene list
pymn.MetaNeighborUS(adata, study_col, bin_col, fast_version=slow)

if not Path(output_path).is_dir():
    logging.info(f'Creating Path {output_path}')
    Path(output_path).mkdir(parents=True)

if rank_by is not None:
    dataset_name = dataset_name + '_ranked'
if join_bin is not None:
    dataset_name = dataset_name + f'_joined_{join_bin}'

term = term.replace(':', '_')

out_fn = f'{output_path}MNUS_{term}_{dataset_name}_{n_bins}.csv'

logging.info(f'Saving Results in: \n {out_fn}')

output_mn_bin_lineage_res(adata.uns['MetaNeighborUS'], out_fn)
