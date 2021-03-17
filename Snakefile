import os
import numpy as np
import pandas as pd

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
## Snakefile for computing metaneighbor and continuity scores
def get_genelist_fn(wildcards):
	fn = config['prepare_terms']['terms_dir'] + wildcards.term.replace(' ','_') + '_genelist.txt'
	print(fn)
	return fn


def build_save_dataset_name(dataset_name, rank_by, join_bin):
	sdn = dataset_name
	if rank_by is not None:
		sdn+= '_ranked'
	if join_bin is not None:
		sdn += f'_joined_{join_bin}'
	return sdn




homedir = os.path.expanduser('~/')
logging.info(homedir)
configfile: 'config.yml'

dataDir = homedir + config['dataDir']

dataset_name = config['dataset_name']
n_bins=config['n_bins']
save_dataset_name = build_save_dataset_name(dataset_name, 
	config['metaneighbor']['rank_by'], 
	config['metaneighbor']['join_bin'])
suffix = f'{save_dataset_name}_{n_bins}.csv'

terms = pd.read_csv(config['prepare_terms']['terms_design'], index_col=0)

term_s = np.sum(terms.values,axis=0)
size_slice = (term_s >= config['prepare_terms']['min_size']) & (term_s <= config['prepare_terms']['max_size'])
terms = terms.columns[size_slice]
assert terms.shape[0] > 0, "Terms must be greater than 0"
fix_terms = np.vectorize(lambda x: x.replace(' ','_'))
terms = fix_terms(terms)
terms_dir =config['prepare_terms']['terms_dir']

terms = terms[:1]
logging.info(terms)
pvers = '/home/bharris/miniconda3/envs/python3_base/bin/python3 -W ignore '

rule all:
	input:
		expand('{dataDir}{dataset}/MNUS_{term}_{suffix}.{extension}', 
			dataset = config['dataset_name'], 
			term=terms, 
			extension=['csv','pdf'],
			suffix=suffix[:-4],
			dataDir=dataDir),
		expand('{dataDir}{dataset}/CONTINUITY_{term}_{suffix}', 
			dataset = config['dataset_name'], 
			term=terms, 
			suffix=suffix,
			dataDir=dataDir),
		expand(terms_dir + '{term}_genelist.txt', term=terms)

rule prepare_terms:
    input:
        term_design=config['prepare_terms']['terms_design']
    params:
        terms_dir=terms_dir,
        min_size=config['prepare_terms']['min_size'],
        max_size=config['prepare_terms']['max_size'],
        terms=fix_terms(terms)
    output:
        expand(terms_dir + '{term}_genelist.txt',term=terms)
    conda:
    	'envs/pseudotime_metaneighbor.yml'
    script:
    	'scripts/prepare_gene_lists.py'

rule metaneighbor:
	input:
		adata_fn=config['dataset_filename'],
		term_fn=get_genelist_fn
	params:
		join_bin = config['metaneighbor']['join_bin'],
		rank_by = config['metaneighbor']['rank_by'],
		dataset_name = config['dataset_name'],
		n_bins=config['n_bins'],
		dataDir=homedir + config['dataDir'],
		study_col=config['metaneighbor']['study_col'],
		pseudotime_col=config['metaneighbor']['pseudotime_col'],
		max_genes=config['metaneighbor']['max_genes'],
		fast_version=config['metaneighbor']['fast_version']
	threads: config['mkl_threads']
	output:
		dataDir + dataset_name + '/MNUS_{term}_' + suffix
	conda:
		'envs/pseudotime_metaneighbor.yml'
	script:
		'scripts/test_term_on_pseudotime.py'
		

rule continuity:
	input:
		mn_res=dataDir + dataset_name + '/MNUS_{term}_' + suffix
	params:
		threshold=config['continuity']['threshold']
	output:
		dataDir + dataset_name + '/CONTINUITY_{term}_' + suffix
	conda:
		'envs/pseudotime_metaneighbor.yml'
	script:
		'scripts/continunity_test_on_term.py'

rule heatmap:
	input:
		mn_res=dataDir + dataset_name + '/MNUS_{term}_' + suffix
	output:
		dataDir + dataset_name + '/MNUS_{term}_' + suffix[:-3] + 'pdf'
	conda:
		'envs/pseudotime_metaneighbor_r.yml'
	# log:
	# 	'log/heatmap_{term}.log'
	script:
		'scripts/metaneighbor_pseudotime_heatmap.R'
	# shell:
	# 	'Rscript scripts/metaneighbor_pseudotime_heatmap.R --fn {input.mn_res}'