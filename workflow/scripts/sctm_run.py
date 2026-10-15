import os
import sys
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import numpy as np
np.random.seed(88)
import scanpy as sc
import squidpy as sq
import torch

# Point rpy2 at the R installation shipped inside the conda env so mclust is
# importable.
os.environ.setdefault('R_HOME', os.path.join(sys.prefix, 'lib', 'R'))

import sctm

N_TOPICS = 20


def mclust_R(embedding, num_cluster, modelNames='EEE', random_seed=2023):
    """Cluster an embedding with the R mclust package via rpy2."""
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library('mclust')
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    robjects.r['set.seed'](random_seed)
    rmclust = robjects.r['Mclust']
    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(embedding),
                  num_cluster, modelNames)
    return np.array(res[-2]).astype(int)


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'scTM', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()
    adata.var_names_make_unique()

    ################################ run method ###############################
    print('\nRunning method')

    print('Preprocessing')
    sc.pp.filter_genes(adata, min_cells=3)
    try:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=3000)
    except ValueError:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3',
                                    n_top_genes=3000, span=0.5)

    print('Building spatial graph')
    sq.gr.spatial_neighbors(adata,
                            coord_type='generic',
                            n_neighs=tech.n_neighbours)

    print('Training STAMP')
    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    model = sctm.stamp.STAMP(adata[:, adata.var.highly_variable],
                             n_topics=N_TOPICS)
    model.train(device=device)

    # STAMP yields an interpretable topic embedding but no domain labels;
    # cluster the cell-by-topic matrix to obtain spatial domains, as the other
    # embedding-based methods in this benchmark do.
    topic_prop = model.get_cell_by_topic()

    print('Clustering')
    labels = mclust_R(topic_prop.to_numpy(), num_cluster=n_domains)

    ############################### save output ###############################

    assert len(labels) == len(coords)

    import pandas as pd
    pd.Series(labels, name='region').to_csv(
        os.path.join(output_folder, 'clustering.tsv'),
        sep='\t',
        index=False)
    print('\n Saved domain labels to clustering.tsv')