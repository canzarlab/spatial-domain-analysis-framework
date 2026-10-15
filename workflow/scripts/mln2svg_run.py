import os
import sys
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import numpy as np
np.random.seed(88)
import random
random.seed(88)
import scanpy as sc
import torch
torch.manual_seed(88)
torch.cuda.manual_seed(88)

# Point rpy2 at the R installation shipped inside the conda env so mclust is
# importable (the upstream config.py hardcodes author-specific paths).
os.environ.setdefault('R_HOME', os.path.join(sys.prefix, 'lib', 'R'))

# MLN2SVG uses flat, relative imports (e.g. `from utils import ...`), so the
# package directory has to be on sys.path rather than imported as a package.
sys.path.insert(0,
                os.path.join('workflow', 'scripts', 'MLN2SVG', 'MLN2SVG'))
from MLN2SVG import MLN2SVG


def preprocess(adata):
    try:
        sc.pp.highly_variable_genes(adata,
                                    flavor='seurat_v3',
                                    n_top_genes=3000)
    except ValueError:
        sc.pp.highly_variable_genes(adata,
                                    flavor='seurat_v3',
                                    n_top_genes=3000,
                                    span=0.5)
    return adata


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'MLN2SVG', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    print('Preprocessing')
    adata = preprocess(adata)

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    mln = MLN2SVG(adata,
                  num_cluster=n_domains,
                  device=device)

    print('Reconstructing spatial expression')
    adata_processed, exp_original, exp_augmented, graph = mln.load_data()
    print('Training')
    latent_original, _, _, _ = mln.train_MLN2SVG(exp_original,
                                                 exp_augmented,
                                                 graph)
    print('Clustering')
    adata_processed = mln.compute_clusters(adata_processed, latent_original)

    ############################### save output ###############################

    labels = adata_processed.obs['mln2svg_clust']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
