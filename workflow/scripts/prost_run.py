import os
import sys
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import scanpy as sc
import torch

# Point rpy2 at the R installation shipped inside the conda env so mclust is
# importable (the upstream example hardcodes an env-specific R_HOME).
# os.environ.setdefault('R_HOME', os.path.join(sys.prefix, 'lib', 'R'))

from PROST import PROST

# Set seed
SEED = 818
PROST.setup_seed(SEED)

if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'PROST', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    platform = dataset.split('_')[0]

    coords_span = (coords.abs().max()-coords.abs().min()).mean()
    if coords_span < 2000:
        adata.obsm['spatial'] = adata.obsm['spatial'] * (2500 / coords_span)
        print('scaling coordinates.')

    ################################ run method ###############################
    print('\nRunning method')

    print(adata)
    print('Calculating PROST Index')
    adata = PROST.prepare_for_PI(adata, platform=platform)
    print(adata)
    adata = PROST.cal_PI(adata, platform=platform)

    print(adata)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata = PROST.feature_selection(adata, by='prost', n_top_genes=3000)

    print('Running PROST Neural Network')
    PROST.run_PNN(adata,
                  key_added='PROST',
                  init='mclust',
                  n_clusters=n_domains,
                  adj_mode='neighbour',
                  k_neighbors=7,
                  lap_filter=2,
                  lr=0.1,
                  SEED=SEED,
                  max_epochs=500,
                  tol=5e-3,
                  post_processing=True,
                  pp_run_times=3,
                  cuda=torch.cuda.is_available())

    ############################### save output ###############################

    labels = adata.obs['pp_clustering']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
