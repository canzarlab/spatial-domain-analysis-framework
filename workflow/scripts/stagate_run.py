import os
from read_input import read_input
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import numpy as np
np.random.seed(88)
import scanpy as sc
from sklearn.neighbors import kneighbors_graph
import scipy.sparse as sp

import tensorflow as tf
tf.compat.v1.disable_eager_execution()
tf.random.set_seed(88)

import STAGATE


def mclust_R(adata, num_cluster, modelNames='EEE', random_seed=2020):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """

    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.X), num_cluster, modelNames)
    try:
        mclust_res = np.array(res[-2])
    except Exception as e:
        print(e)
        mclust_res = [0] * adata.X.shape[0]

    adata.obs['mclust'] = mclust_res
    adata.obs['mclust'] = adata.obs['mclust'].astype('int')
    adata.obs['mclust'] = adata.obs['mclust'].astype('category')
    return adata


def preprocess(adata):
    try:
        sc.pp.highly_variable_genes(adata,
                                    flavor="seurat_v3",
                                    n_top_genes=3000)
    except ValueError:
        sc.pp.highly_variable_genes(adata,
                                    flavor="seurat_v3",
                                    n_top_genes=3000,
                                    span=0.5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata


def graph(adata, k_neighbors):
    adata.obsm['adj'] = kneighbors_graph(adata.obsm['spatial'],
                                         n_neighbors=k_neighbors)
    return adata.obsm['adj']


def train(adata, adj):
    adj_sp = sp.csr_matrix(adj)

    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']

    KNN_list = []
    for rowi in range(len(adj_sp.indptr)-1):
        indices = adj_sp.indices[adj_sp.indptr[rowi]:adj_sp.indptr[rowi+1]]
        data = adj_sp.data[adj_sp.indptr[rowi]:adj_sp.indptr[rowi+1]]
        KNN_list.append(pd.DataFrame(zip([rowi] * len(indices),
                                         indices,
                                         data)))

    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']

    Spatial_Net = KNN_df.copy()
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance'] > 0,]
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index), ))
    Spatial_Net['Cell1'] = Spatial_Net['Cell1'].map(id_cell_trans)
    Spatial_Net['Cell2'] = Spatial_Net['Cell2'].map(id_cell_trans)

    adata.uns['Spatial_Net'] = Spatial_Net
    adata = STAGATE.train_STAGATE(adata)
    adata.obsm['embedding'] = adata.obsm['STAGATE']

    return adata.obsm['embedding']


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'STAGATE', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    print('Preprocessing')
    adata = preprocess(adata)
    print('Calculating graph')
    adj = graph(adata, k_neighbors=tech.n_neighbours)
    print('Training')
    _ = train(adata, adj=adj)

    adata = STAGATE.mclust_R(adata, num_cluster=n_domains)

    ############################### save output ###############################

    labels = adata.obs['mclust']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
