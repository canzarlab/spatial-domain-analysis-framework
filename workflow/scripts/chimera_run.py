#### SEDR

from sklearn.cluster import KMeans
import os
import argparse
from tqdm import tqdm
import subprocess
import pandas as pd
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

def SEDR_preprocess(adata):
    print('\nPreprocessing\n')
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
    adata_new = sc.AnnData(adata_X)
    adata_new.obsm['spatial'] = adata.obsm['spatial'].copy()
    adata_new.var['highly_variable'] = [True] * adata_new.X.shape[1]

    return adata_new

def SEDR_adjacency(adata):
    print('\nRunning method')
    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    print('===== Using device: ' + device)
    
    SEDR.fix_seed(2023)
    '''
    graph_dict = SEDR.graph_construction(adata)
    '''
    return sp.csr_matrix(SEDR.graph_func.generate_adj_mat(adata, include_self=False, n=6))

def SEDR_train(adata, adj):
    adj_m1 = adj
    adj_m1 = sp.coo_matrix(adj_m1)

    # Store original adjacency matrix (without diagonal entries) for later
    adj_m1 = adj_m1 - sp.dia_matrix((adj_m1.diagonal()[np.newaxis, :], [0]), shape=adj_m1.shape)
    adj_m1.eliminate_zeros()

    # Some preprocessing
    adj_norm_m1 = SEDR.graph_func.preprocess_graph(adj_m1)
    adj_m1 = adj_m1 + sp.eye(adj_m1.shape[0])
    # adj_label_m1 = torch.FloatTensor(adj_label_m1.toarray())

    adj_m1 = adj_m1.tocoo()
    shape = adj_m1.shape
    values = adj_m1.data
    indices = np.stack([adj_m1.row, adj_m1.col])
    adj_label_m1 = torch.sparse_coo_tensor(indices, values, shape)

    norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - adj_m1.sum()) * 2)

    # # generate random mask
    # adj_mask = mask_generator(adj_label_m1.to_sparse(), N)

    graph_dict = {
        "adj_norm": adj_norm_m1,
        "adj_label": adj_label_m1.coalesce(),
        "norm_value": norm_m1,
        # "mask": adj_mask
    }
    sedr_net = SEDR.Sedr(adata.X,
                         graph_dict,
                         device=device)
    sedr_net.train_with_dec()

    sedr_feat, _, _, _ = sedr_net.process()
    adata.obsm['embedding'] = sedr_feat

    return adata

#### CCST

import os
import argparse
from tqdm import tqdm
import subprocess
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
    pca_n_comps = min(adata.n_vars, adata.n_obs, 200)
    adata_X = sc.pp.pca(adata_X,
                        n_comps=pca_n_comps)
    return adata_X


def CCST_preprocess(adata):
    features = adata_preprocess(adata,
                                min_cells=3)
    adata_new = sc.AnnData(features)
    adata_new.obsm['spatial'] = adata.obsm['spatial'].copy()
    adata_new.var['highly_variable'] = [True] * adata_new.X.shape[1]
    return adata_new


def CCST_adjacency(adata, k_lo=2.8, k_hi=9.8, lambda_I=.3, lo=0, hi=100000):
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


def CCST_train(adata, adj, lambda_I, output_folder, device=('cuda' if torch.cuda.is_available() else 'cpu')):
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
    return adata


#### SCANIT

import os
import argparse
from tqdm import tqdm
import subprocess
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

import scanit
import scanpy as sc
from somde import SomNode
import numpy as np
import networkx as nx
import gudhi
import torch


def SCANIT_preprocess(adata):
    n_sv_genes = 3000
    adata_sp = adata
    sc.pp.normalize_total(adata_sp)
    pts = adata.obsm['spatial']
    df_sp = pd.DataFrame(data=adata.X.copy(), columns=list(adata.var_names))
    som = SomNode(pts, 5)
    ndf, ninfo = som.mtx(df_sp.T)
    _ = som.norm()
    result, SVnum = som.run()

    df_somde = result
    sv_genes = list(df_somde['g'].values[:n_sv_genes])
    adata_sp = adata_sp[:, sv_genes]
    sc.pp.log1p(adata_sp)
    sc.pp.scale(adata_sp)
    adata_sp.var['highly_variable'] = [True] * adata_sp.X.shape[1]
    
    return adata_sp


def SCANIT_adjacency(adata, k_neighbors=5):
    pts = adata.obsm['spatial'].copy()
    cut = scanit._utils.estimate_cutoff_knn(pts)
    A = scanit._utils.graph_alpha(pts, cut=cut)

    adata.obsp['scanit-graph'] = A
    adata.obsm['adj'] = adata.obsp['scanit-graph']

    return adata.obsm['adj']


def SCANIT_train(adata, adj, device):
    adata.obsp['scanit-graph'] = sp.csr_matrix(adj)
    adata_sp = adata.copy()
    #adata_sp.X = sp.csr_matrix(adata_sp.X)
    scanit.tl.spatial_representation(adata_sp,
                                     device=device,
                                     python_seed=0,
                                     torch_seed=0,
                                     numpy_seed=0)

    adata.obsm['embedding'] = adata_sp.obsm['X_scanit'].copy()
    return adata


#### STAGATE

import os
import argparse
from tqdm import tqdm
import subprocess
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

def STAGATE_preprocess(adata):
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


def STAGATE_adjacency(adata, k_neighbors=6):
    adata.obsm['adj'] = kneighbors_graph(adata.obsm['spatial'],
                                         n_neighbors=k_neighbors)
    return adata.obsm['adj']


def STAGATE_train(adata, adj):
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

    adata_new = sc.AnnData(sp.csr_matrix(adata.X))
    adata_new.uns['Spatial_Net'] = Spatial_Net
    adata_new.obsm['spatial'] = adata.obsm['spatial'].copy()
    adata_new.var['highly_variable'] = [True] * adata_new.X.shape[1]
    adata = STAGATE.train_STAGATE(adata_new)
    adata.obsm['embedding'] = adata.obsm['STAGATE']

    return adata


### GraphST

from GraphST.preprocess import preprocess as Gpreprocess
from GraphST.preprocess import construct_interaction_KNN
from GraphST.GraphST import GraphST

def GraphST_preprocess(adata):
    Gpreprocess(adata)
    return adata

def GraphST_adjacency(adata):
    construct_interaction_KNN(adata)
    return adata.obsm['adj']

def GraphST_train(adata, adj, datatype='Slide'):
    adata.obsm['adj'] = adj
    adata.obsm['graph_neigh'] = 1*(adj > 0)
    adata.X = sp.csr_matrix(adata.X)
    #print(adata)
    #datatype = 'Slide' if 'Slide' in tech.name else '10x'
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model = GraphST(adata,
                    device=device,
                    datatype=datatype)
    adata = model.train()
    pca = PCA(n_components=20, random_state=42) 
    embedding = pca.fit_transform(adata.obsm['emb'].copy())
    adata.obsm['embedding'] = embedding

    return adata
    
#### SpaceFlow

from SpaceFlow import SpaceFlow

def SpaceFlow_preprocess(adata):
    sc.pp.filter_genes(adata, min_counts=1)
    sf = SpaceFlow.SpaceFlow(adata=adata)
    sf.preprocessing_data(n_top_genes=3000)
    return sf.adata_preprocessed

def SpaceFlow_adjacency(adata):
    adj = SpaceFlow.SpaceFlow.graph_alpha(None, adata.obsm['spatial'])
    adata.obsm['adj'] = adj
    return adata.obsm['adj']

def SpaceFlow_train(adata, adj):
    sf = SpaceFlow.SpaceFlow(adata)
    sf.adata_preprocessed = adata.copy()
    sf.spatial_graph = sp.csr_matrix(adj)
    sf.train(os.path.join(output_folder, 'embedding.txt'))
    adata.obsm['embedding'] = sf.embedding
    return adata
    
#### MAIN()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', '-d', type=str, required=True)
    parser.add_argument('--sample', '-s', type=str, required=True)
    parser.add_argument('--n_domains', '-n', type=int, required=True)
    parser.add_argument('--technology', '-t', type=str)
    meths = [ 'SEDR', 'STAGATE', 'CCST', 'SCANIT', 'GraphST', 'SpaceFlow' ]
    parser.add_argument('--name', type=str, required=True)
    parser.add_argument('--preprocess', type=str, required=True, choices=meths)
    parser.add_argument('--adjacency', type=str, required=True, choices=meths)
    parser.add_argument('--train', type=str, required=True, choices=meths)
    parser.add_argument('--clust', type=str, required=True, choices=['kmeans', 'mclust'])
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
    output_folder = os.path.join('results', dataset, args.name, sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    file = os.path.join(input_folder, f'{sample}_UMImatrix.csv')
    print(file)
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

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    if args.preprocess == 'SEDR': adata = SEDR_preprocess(adata)
    elif args.preprocess == 'CCST': adata = CCST_preprocess(adata)
    elif args.preprocess == 'STAGATE': adata = STAGATE_preprocess(adata)
    elif args.preprocess == 'GraphST': adata = GraphST_preprocess(adata)
    elif args.preprocess == 'SpaceFlow': adata = SpaceFlow_preprocess(adata)
    else: adata = SCANIT_preprocess(adata)

    if args.adjacency == 'SEDR': adj = SEDR_adjacency(adata)
    elif args.adjacency == 'CCST': adj = CCST_adjacency(adata, hi=adata.obsm['spatial'].max())
    elif args.adjacency == 'STAGATE': adj = STAGATE_adjacency(adata)
    elif args.adjacency == 'GraphST': adj = GraphST_adjacency(adata)
    elif args.adjacency == 'SpaceFlow': adj = SpaceFlow_adjacency(adata)
    else: adj = SCANIT_adjacency(adata)

    if args.train == 'SEDR': pred = SEDR_train(adata, adj)
    elif args.train == 'CCST': pred = CCST_train(adata, adj, 0.3, output_folder)
    elif args.train == 'STAGATE': pred = STAGATE_train(adata, adj)
    elif args.train == 'GraphST': pred = GraphST_train(adata, adj)
    elif args.train == 'SpaceFlow': pred = SpaceFlow_train(adata, adj)
    else: pred = SCANIT_train(adata, adj, device=device)


    if args.clust == 'kmeans': pred.obs['pred'] = KMeans(n_clusters=args.n_domains, random_state=42).fit_predict(pred.obsm['embedding'])
    else: mclust_R(pred, args.n_domains, use_rep='embedding', key_added='pred')

    ###### SAVE CLUSTERING
    labels = pred.obs['pred']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
    ### SEDR
