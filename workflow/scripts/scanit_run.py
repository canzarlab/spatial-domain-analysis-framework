import os
from read_input import read_input
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

import scanit
import scanpy as sc
from somde import SomNode
import torch


def preprocess(adata):
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

    return adata_sp


def graph(adata, k_neighbors=5):
    pts = adata.obsm['spatial'].copy()
    cut = scanit._utils.estimate_cutoff_knn(pts)
    A = scanit._utils.graph_alpha(pts, cut=cut)

    adata.obsp['scanit-graph'] = A
    adata.obsm['adj'] = adata.obsp['scanit-graph']

    return adata.obsm['adj']


def train(adata, adj, device):
    adata.obsp['scanit-graph'] = adj
    adata_sp = adata.copy()

    scanit.tl.spatial_representation(adata_sp,
                                     device=device,
                                     python_seed=0,
                                     torch_seed=0,
                                     numpy_seed=0)

    adata.obsm['embedding'] = adata_sp.obsm['X_scanit'].copy()
    return adata_sp.obsm['X_scanit']


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'SCAN-IT', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts.to_numpy())
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    print('Preprocessing')
    adata = preprocess(adata)
    print('Calculating graph')
    adj = graph(adata)
    print('Training')
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    emb = train(adata, adj=adj, device=device)

    adata_emb = sc.AnnData(emb)
    adata_emb.obsm['spatial'] = coords.to_numpy()

    if tech.seq_or_img == 'seq':
        n_neighs = 15
    elif tech.seq_or_img == 'img':
        n_neighs = 10
    else:
        print("no good :(")

    sc.pp.neighbors(adata_emb, n_neighbors=n_neighs)
    lo, hi = 0., 5.
    while hi - lo > 1e-4:
        mi = (lo+hi)/2
        sc.tl.leiden(adata_emb, resolution=mi)
        lein = len(adata_emb.obs['leiden'].unique())
        if lein == n_domains:
            break
        if lein < n_domains:
            lo = mi
        else:
            hi = mi

    ############################### save output ###############################

    labels = adata_emb.obs[['leiden']]
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
