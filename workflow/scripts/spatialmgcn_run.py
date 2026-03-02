from __future__ import division
from __future__ import print_function

import os
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import numpy as np
np.random.seed(88)
import scanpy as sc

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'Spatial-MGCN/Spatial-MGCN/'))

from utils import features_construct_graph, spatial_construct_graph1
import os
from config import Config

import torch.optim as optim
from utils import *
from models import Spatial_MGCN

from read_input import read_input


def train():
    model.train()
    optimizer.zero_grad()
    com1, com2, emb, pi, disp, mean = model(features, sadj, fadj)
    zinb_loss = ZINB(pi, theta=disp, ridge_lambda=0).loss(features, mean, mean=True)
    reg_loss = regularization_loss(emb, graph_nei, graph_neg)
    con_loss = consistency_loss(com1, com2)
    total_loss = config.alpha * zinb_loss + config.beta * con_loss + config.gamma * reg_loss
    emb = pd.DataFrame(emb.cpu().detach().numpy()).fillna(0).values
    mean = pd.DataFrame(mean.cpu().detach().numpy()).fillna(0).values
    total_loss.backward()
    optimizer.step()
    return emb, mean, zinb_loss, reg_loss, con_loss, total_loss

def normalize(adata, highly_genes=3000):
    print("start select HVGs")
    sc.pp.filter_genes(adata, min_cells=100)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=highly_genes)
    adata = adata[:, adata.var['highly_variable']].copy()
    adata.X = adata.X / np.sum(adata.X, axis=1).reshape(-1, 1) * 10000
    sc.pp.scale(adata, zero_center=False, max_value=10)
    return adata



if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'SpatialMGCN', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()


    #### preprocess

    config_file = os.path.join(os.path.dirname(__file__), 'Spatial-MGCN/Spatial-MGCN/config/DLPFC.ini')
    config = Config(config_file)
    config.fdim = min(config.fdim, adata.n_vars)
    highly_genes, k, radius = config.fdim, config.k, config.radius
    
    adata = normalize(adata, highly_genes=highly_genes)
    fadj = features_construct_graph(adata.X, k=k)
    sadj, graph_nei, graph_neg = spatial_construct_graph1(adata, radius=radius)

    adata.obsm["fadj"] = fadj
    adata.obsm["sadj"] = sadj
    adata.obsm["graph_nei"] = graph_nei.numpy()
    adata.obsm["graph_neg"] = graph_neg.numpy()

    features = torch.FloatTensor(adata.X)
    #labels = adata.obs['ground']
    fadj = adata.obsm['fadj']
    sadj = adata.obsm['sadj']
    nfadj = normalize_sparse_matrix(fadj + sp.eye(fadj.shape[0]))
    nfadj = sparse_mx_to_torch_sparse_tensor(nfadj)
    nsadj = normalize_sparse_matrix(sadj + sp.eye(sadj.shape[0]))
    nsadj = sparse_mx_to_torch_sparse_tensor(nsadj)
    graph_nei = torch.LongTensor(adata.obsm['graph_nei'])
    graph_neg = torch.LongTensor(adata.obsm['graph_neg'])

    sadj = nsadj ##
    fadj = nfadj ##

    cuda = not config.no_cuda and torch.cuda.is_available()
    use_seed = not config.no_seed
    
    #_, ground = np.unique(np.array(labels, dtype=str), return_inverse=True)
    #ground = torch.LongTensor(ground)
    config.n = len(counts)#len(ground)
    config.class_num = n_domains#len(ground.unique())

    config.epochs = 200
    config.epochs = config.epochs + 1

    if cuda:
        features = features.cuda()
        sadj = sadj.cuda()
        fadj = fadj.cuda()
        graph_nei = graph_nei.cuda()
        graph_neg = graph_neg.cuda()



    model = Spatial_MGCN(nfeat=config.fdim,
                             nhid1=config.nhid1,
                             nhid2=config.nhid2,
                             dropout=config.dropout)
    if cuda:
        model.cuda()
    optimizer = optim.Adam(model.parameters(), lr=config.lr, weight_decay=config.weight_decay)

    epoch_max = 0
    ari_max = 0
    idx_max = []
    mean_max = []
    emb_max = []

    
    for epoch in range(config.epochs):
        emb, mean, zinb_loss, reg_loss, con_loss, total_loss = train()
        print(dataset, ' epoch: ', epoch, ' zinb_loss = {:.2f}'.format(zinb_loss),
              ' reg_loss = {:.2f}'.format(reg_loss), ' con_loss = {:.2f}'.format(con_loss),
              ' total_loss = {:.2f}'.format(total_loss))
        kmeans = KMeans(n_clusters=config.class_num).fit(emb)
        idx = kmeans.labels_

    adata.obsm['emb'] = emb
    adata.obsm['mean'] = mean
    adata.obs['labels'] = idx
    
    labels = adata.obs[['labels']]
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')

