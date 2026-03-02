import os
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import scanpy as sc
from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable. (<- comment from SEDR tutorial)
import torch
import rpy2

import SEDR


def mclust_R(adata, n_domains, use_rep='SEDR', key_added='SEDR', random_seed=2023):  # function taken directly from SEDR, minus indicated line
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    import os
    # os.environ['R_HOME'] = '/scbio4/tools/R/R-4.0.3_openblas/R-4.0.3'  # <- this line removed wrt SEDR
    modelNames = 'EEE'

    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm[use_rep]), n_domains, modelNames)
    mclust_res = np.array(res[-2])

    adata.obs[key_added] = mclust_res
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    return adata


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'SEDR', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ############################# preprocess data #############################
    print('\nPreprocessing\n')

    # preprocessing taken from SEDR.adata_preprocess defaults.
    # function can't be used due to AttributeError for adata.X.toarray()
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.filter_genes(adata, min_counts=10)
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.highly_variable_genes(adata,
                                flavor="seurat_v3",
                                n_top_genes=2000,
                                subset=True)
    sc.pp.scale(adata)

    n_comp = min(adata.n_vars, adata.n_obs, 200)
    adata_X = PCA(n_components=n_comp,
                  random_state=42
                  ).fit_transform(adata.X)
    adata.obsm['X_pca'] = adata_X

    ################################ run method ###############################
    print('\nRunning method')

    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    print('===== Using device: ' + device)

    SEDR.fix_seed(2023)
    graph_dict = SEDR.graph_construction(adata)
    sedr_net = SEDR.Sedr(adata.obsm['X_pca'],
                         graph_dict,
                         device=device)
    sedr_net.train_with_dec()

    sedr_feat, _, _, _ = sedr_net.process()
    adata.obsm['SEDR'] = sedr_feat

    mclust_R(adata, n_domains)

    ############################### save output ###############################

    labels = adata.obs['SEDR']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
