import os
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import torch
torch.manual_seed(88)
torch.cuda.manual_seed(88)
import scanpy as sc

import mnmstpy as mnmst


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'MNMST', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=3000)
    hvg_filter = adata.var['highly_variable']
    sc.pp.normalize_total(adata, inplace=True)
    adata = adata[:, hvg_filter]
    enhanced_adata, cell_spatial = mnmst.data_enhance(adata, k_nei=6)
    sc.pp.pca(enhanced_adata)
    low_dim_x = enhanced_adata.obsm['X_pca']

    init_graph = mnmst.construct_initial_graph(low_dim_x)
    
    device = torch.device('cuda') if torch.cuda.is_available() else 'cpu'
    init_W_tensor = torch.from_numpy(init_graph).float().to(device)
    cell_spatia_tensor = torch.from_numpy(cell_spatial.toarray()).float().to(device)
    x_tensor = torch.from_numpy(low_dim_x.T).to(device)
    C_gpu = mnmst.sparse_self_representation_torch(x_tensor, init_w=init_W_tensor, device=device)

    # datatype is left at default "10x" for all except Slide-seq data
    datatype = 'Slide' if 'Slide' in tech.name else '10x'
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    Z_gpu = mnmst.MNMST_representation_gpu(C_gpu, cell_spatia_tensor, device=device)
    Z = Z_gpu.detach().cpu().numpy()

    key_added = 'representation'
    conns_key = 'representation'
    dists_key = 'representation'

    adata.uns[key_added] = {}    
    representation_dict = adata.uns[key_added]
    representation_dict['connectivities_key'] = conns_key
    representation_dict['distances_key'] = dists_key
    representation_dict['var_names_use'] = adata.var_names.to_numpy()
    
    representation_dict['params'] = {}
    representation_dict['params']['method'] = 'umap'
    adata.obsp['representation'] = Z

    lo, hi = 0., 5.
    while hi - lo > 1e-4:
        mi = (lo+hi)/2
        sc.tl.leiden(adata, neighbors_key='representation', resolution=mi, key_added='MNMST_GPU')
        lein = len(adata.obs['MNMST_GPU'].unique())
        if lein == n_domains:
            break
        if lein < n_domains:
            lo = mi
        else:
            hi = mi

    ############################### save output ###############################

    labels = adata.obs['MNMST_GPU']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
