import os
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import scanpy as sc
import SpaGCN as spg
import random
import torch


def preprocess(adata):
    spg.prefilter_genes(adata)
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata


def graph(adata):
    adj = spg.calculate_adj_matrix(x=adata.obs['x'],
                                   y=adata.obs['y'],
                                   histology=False)

    adata.obsm['adj'] = adj
    return adj


def train(adata, n_domains, adj):
    r_seed = t_seed = n_seed = 66

    l = spg.search_l(p=0.5, adj=adj)
    assert l is not None

    res = spg.search_res(adata=adata,
                         adj=adj,
                         l=l,
                         target_num=n_domains)

    clf = spg.SpaGCN()
    clf.set_l(l)

    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)

    clf.train(adata=adata, adj=adj, res=res)
    adata.obs["pred"], _ = clf.predict()

    emb = clf.model.predict(clf.embed, clf.adj_exp)[0].detach().cpu()
    return emb.numpy()


def refine(adata):
    adj_2d = spg.calculate_adj_matrix(x=adata.obs['x'],
                                      y=adata.obs['y'],
                                      histology=False)
    refined_pred = spg.refine(sample_id=adata.obs.index.tolist(),
                              pred=adata.obs['pred'].tolist(),
                              dis=adj_2d)
    adata.obs["refined_pred"] = refined_pred

    return adata


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'SpaGCN', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts, obs=coords)

    ################################ run method ###############################
    print('\nRunning method')

    print('Preprocessing')
    adata = preprocess(adata)
    print('Calculating graph')
    adj = graph(adata)
    print('Training')
    _ = train(adata, n_domains, adj=adj)
    print('Refining')
    adata = refine(adata)

    ############################### save output ###############################

    labels = adata.obs['refined_pred']
    if len(labels) != len(coords):
        idx_missing = list(set(map(str, coords.index)) - set(labels.index))
        print(f'Missing {len(idx_missing)} cells. Inserting NaN.')
        for idx in idx_missing:
            labels.loc[idx] = np.nan

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
