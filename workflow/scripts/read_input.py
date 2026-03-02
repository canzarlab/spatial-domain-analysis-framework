import argparse
import pandas as pd
import os
from tqdm import tqdm
import subprocess


def read_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', '-d', type=str, required=True)
    parser.add_argument('--sample', '-s', type=str, required=True)
    parser.add_argument('--n_domains', '-n', type=int, required=True)
    parser.add_argument('--technology', '-t', type=str)
    args = parser.parse_args()
    dataset = args.dataset
    sample = args.sample
    n_domains = args.n_domains
    tech = args.technology if args.technology else dataset.split('_')[0]
    tech = pd.read_csv(os.path.join('config', 'technology_meta.csv'),
                       index_col=0).loc[tech]
    assert tech.seq_or_img in ['seq', 'img']

    ################################ read data ################################
    print('Reading data')

    input_folder = os.path.join('resources', dataset)
    file = os.path.join(input_folder, f'{sample}_UMImatrix.csv')
    n_lines = int(subprocess.run(['wc', '-l', file],
                                 capture_output=True,
                                 text=True
                                 ).stdout.split(' ')[0])
    counts = pd.concat([chunk for chunk
                        in tqdm(pd.read_csv(file, header=None, chunksize=100),
                                desc='count matrix progress',
                                total=round(n_lines/100))])
    coords = pd.read_csv(os.path.join(input_folder, f'{sample}_coords.csv'),
                         header=None,
                         names=['x', 'y'])
    assert counts.shape[0] == coords.shape[0]
    return dataset, sample, n_domains, tech, counts, coords
