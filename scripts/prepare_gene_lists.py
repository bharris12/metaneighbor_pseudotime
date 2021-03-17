import pandas as pd
import numpy as np
from pathlib import Path

from argparse import ArgumentParser

import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

from types import SimpleNamespace

if 'snakemake' in globals():
    logging.info("READING PARAMS FROM SNAKEMAKE")
    args = dict(design_fn=snakemake.input[0],
                output_dir=snakemake.params[0],
                min_size=snakemake.params[1],
                max_size=snakemake.params[2])
    args = SimpleNamespace(**args)
else:
    logging.info("READING PARAMS FROM CMD LINE")
    parser = ArgumentParser(description='Get Stuff')

    parser.add_argument(
        '--design-fn',
        type=str,
        required=True,
        help='File name for design matrix, must be genes x terms')
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help=
        'Directory to output the saved gene lists and terms x filenames file to'
    )
    parser.add_argument('--min-size',
                        type=int,
                        default=10,
                        help='Minimum size of gene list')
    parser.add_argument('--max-size',
                        type=int,
                        default=10,
                        help='Maximum size of gene list')

    args = parser.parse_args()

design = pd.read_csv(args.design_fn, index_col=0)

output_path = args.output_dir
if not Path(output_path).is_dir():
    logging.info(f'Creating Path {output_path}')
    Path(output_path).mkdir()

term_sizes = design.sum()
design = design.loc[:, (term_sizes > args.min_size) &
                    (term_sizes < args.max_size)]

assert design.shape[
    1] > 0, f'No design are larger than {args.min_size} and smaller than {args.max_size}'

logging.info(f'Creating gene list files for {design.shape[1]} terms')

term_files = pd.Series(index=design.columns, name='terms_file')
term_files.index.name = 'term_name'

for term in design.columns:
    genes = design.index[design[term].astype(bool)]
    term = term.replace(' ', '_')
    outfile = f'{output_path}{term}_genelist.txt'
    np.savetxt(outfile, genes, fmt='%s')
    term_files[term] = outfile

term_files.to_csv(f'{output_path}terms_info.csv')
