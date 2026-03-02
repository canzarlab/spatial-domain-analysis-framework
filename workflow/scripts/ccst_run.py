import os
from read_input import read_input
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
import numpy as np
from scipy import sparse
import torch
torch.manual_seed(88)
torch.cuda.manual_seed(88)
from torch_geometric.data import DataLoader

from CCST.CCST import get_graph, train_DGI, train_DGI, Kmeans_cluster


def adata_preprocess(i_adata, min_cells=3):
    sc.pp.filter_genes(i_adata,
                       min_cells=min_cells)
    adata_X = sc.pp.normalize_total(i_adata,
                                    target_sum=1,
                                    exclude_highly_expressed=True,
                                    inplace=False)['X']
    adata_X = sc.pp.scale(adata_X)
    pca_n_comps = min(adata.n_vars-1, adata.n_obs-1, 200)
    adata_X = sc.pp.pca(adata_X,
                        n_comps=pca_n_comps)
    return adata_X


def preprocess(adata):
    features = adata_preprocess(adata,
                                min_cells=3)
    adata_new = sc.AnnData(features)
    adata_new.obsm['spatial'] = adata.obsm['spatial'].copy()
    adata_new.var['highly_variable'] = [True] * adata_new.X.shape[1]
    return adata_new


def graph(adata, k_lo=2.8, k_hi=9.8, lambda_I=.3, lo=0, hi=100000):
    '''
    k_lo, k_hi is the range in which we want our avg neighbour count to be in
    lo, hi is the starting range for bisection
    '''
    coordinates = adata.obsm['spatial']
    cell_num = len(coordinates)

    print('calculating distance matrix, it takes a while')
    tcord = torch.from_numpy(coordinates)
    distance_array = torch.cdist(tcord, tcord).numpy()

    threshold = (lo+hi)/2
    while hi - lo > 1e-3:
        threshold = (lo+hi)/2
        num_big = np.where(distance_array < threshold)[0].shape[0]

        avg_num = num_big/(cell_num*2)
        print('avg num', avg_num)
        if k_lo <= avg_num <= k_hi:
            break
        if avg_num < k_lo:
            lo = threshold
        else:
            hi = threshold

    distance_matrix = distance_array
    distance_matrix_threshold_I = np.zeros(distance_matrix.shape)
    distance_matrix_threshold_W = np.zeros(distance_matrix.shape)
    for i in range(distance_matrix_threshold_I.shape[0]):
        for j in range(distance_matrix_threshold_I.shape[1]):
            if distance_matrix[i, j] <= threshold and distance_matrix[i, j] > 0:
                distance_matrix_threshold_I[i, j] = 1
                distance_matrix_threshold_W[i, j] = distance_matrix[i, j]

    ############### get normalized sparse adjacent matrix
    distance_matrix_threshold_I_N = np.float32(distance_matrix_threshold_I)
    distance_matrix_threshold_I_N_crs = sparse.csr_matrix(distance_matrix_threshold_I_N)

    adj_0 = distance_matrix_threshold_I_N_crs

    num_points = adata.X.shape[0]
    adj_I = np.eye(num_points)
    adj_I = sparse.csr_matrix(adj_I)
    adj = (1-lambda_I)*adj_0 + lambda_I*adj_I

    adata.obsm['adj'] = adj

    return adata.obsm['adj']


class Params:
    def __init__(self, **entries):
        self.__dict__.update(entries)


def train(adata, adj, lambda_I, output_folder, device=('cuda' if torch.cuda.is_available() else 'cpu')):
    args = Params(
        hidden=256,
        load=False,
        num_epoch=5000,
        lambda_I=lambda_I,
        model_path=output_folder
        )

    adj = sparse.csr_matrix(adj)
    data_list = get_graph(adj, adata.X.copy())
    data_loader = DataLoader(data_list, batch_size=1)
    DGI_model = train_DGI(args,
                          data_loader=data_loader,
                          in_channels=adata.n_vars)

    for data in data_loader:
        data.to(device)
        X_embedding, _, _ = DGI_model(data)
        X_embedding = X_embedding.cpu().detach().numpy()

    adata.obsm['embedding'] = X_embedding
    return adata.obsm['embedding']


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'CCST', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy().astype('float')

    ################################ run method ###############################
    print('\nRunning CCST')

    print('Preprocessing')
    adata = preprocess(adata)
    print('Calculating graph')
    if tech.seq_or_img == 'seq':
        lambda_I = 0.3
        k_lo = 2
        k_hi = 5
    elif tech.seq_or_img == 'img':
        lambda_I = 0.8
        k_lo = 13
        k_hi = 17
    else:
        raise Exception('this part should not be reached.')

    adj = graph(adata, k_lo=k_lo, k_hi=k_hi,
                lambda_I=lambda_I,
                hi=coords.max(axis=None))
    print('Training')
    emb = train(adata, adj=adj, lambda_I=lambda_I, output_folder=output_folder)

    cluster_labels, _ = Kmeans_cluster(emb, n_domains)

    ############################### save output ###############################

    labels = pd.DataFrame(cluster_labels)
    assert len(labels) == len(coords)
    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
