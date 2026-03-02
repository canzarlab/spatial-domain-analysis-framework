import os
import argparse
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
os.environ['OMP_NUM_THREADS'] = '2'
os.environ['MKL_NUM_THREADS'] = '2'
os.environ['NUMEXPR_NUM_THREADS'] = '2'
from pathlib import Path
import numpy as np
np.random.seed(0)
from sklearn.preprocessing import StandardScaler
import torch
torch.set_num_threads(4)

from SpiceMix.SpiceMix.model import SpiceMix
from SpiceMix.SpiceMix.util import smooth_discrete, smooth_continuous, clustering_louvain_nclust


# this function is unchanged from SpiceMix example notebook
# (except for generalising n_domains)
def calculate_clustering(latent_state: np.ndarray, df_meta: pd.DataFrame, num_nbrs):
    latent_state = StandardScaler().fit_transform(latent_state)
    latent_state = smooth_continuous(df_meta[['x', 'y']].values, latent_state, df_meta['repli'])
    cluster = clustering_louvain_nclust(
        latent_state.copy(), n_domains,
        kwargs_neighbors=dict(n_neighbors=num_nbrs),
        kwargs_clustering=dict(),
        force_n_cluster=False
    )
    while True:
        cluster_prev = cluster
        cluster = smooth_discrete(df_meta[['x', 'y']].values, cluster, df_meta['repli'])
        if np.all(cluster_prev == cluster):
            break
    return cluster


if __name__ == '__main__':
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

    ################################ run method ###############################
    print('\nRunning method')

    data_folder = Path(os.path.join('results', dataset, 'SpiceMix', sample))
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)

    repli_list = [sample]

    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    # default arguments from main.py
    obj = SpiceMix(K=15,  # default
                   lambda_Sigma_x_inv=1e-6,
                   power_Sigma_x_inv=2,
                   repli_list=[sample],
                   context=dict(device=device, dtype=torch.float64),
                   )

    obj.load_dataset(data_folder)
    obj.initialize()

    for iiter in range(10):  # default from main.py
        obj.estimate_weights(iiter=iiter, use_spatial=[False]*obj.num_repli)
        obj.estimate_parameters(iiter=iiter, use_spatial=[False]*obj.num_repli)

    obj.initialize_Sigma_x_inv()
    for iiter in range(1, 201):   # default from main.py
        obj.estimate_parameters(iiter=iiter, use_spatial=[True]*obj.num_repli)
        obj.estimate_weights(iiter=iiter, use_spatial=[True]*obj.num_repli)

    df = pd.read_csv(os.path.join(data_folder, 'files', f'meta_{sample}.csv'))
    df['repli'] = sample

    latent_state = obj.Xs[0].cpu().numpy()

    cluster = calculate_clustering(latent_state, df, num_nbrs=20)  # value from example notebook

    ############################### save output ###############################

    labels = pd.DataFrame(cluster)
    assert len(labels) == len(df)

    labels.to_csv(os.path.join(data_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
