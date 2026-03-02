import os
from read_input import read_input
import numpy as np
import scanpy as sc
from Banksy_py.banksy_utils.filter_utils import normalize_total, filter_hvg
from Banksy_py.banksy.initialize_banksy import initialize_banksy
from Banksy_py.banksy.embed_banksy import generate_banksy_matrix
from Banksy_py.banksy_utils.umap_pca import pca_umap
from Banksy_py.banksy.cluster_methods import run_mclust_partition
import warnings
warnings.filterwarnings('ignore')


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'BANKSY', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts, obs=coords)

    ############################# preprocess data #############################
    print('\nPreprocessing\n')

    # note: not using BANKSYs filter_cells function because it filters mito-
    # chondrial genes by default, while we're not providing gene information
    print(f'cells before filtering: {adata.n_obs}')
    sc.pp.filter_cells(adata, min_counts=40)
    # sc.pp.filter_cells(adata, max_counts=1000)
    print(f'cells after filtering: {adata.n_obs}')
    print(f'genes before filtering: {adata.n_vars}')
    sc.pp.filter_genes(adata, min_cells=10)
    print(f'genes after filtering: {adata.n_vars}')

    adata = normalize_total(adata)
    adata, adata_allgenes = filter_hvg(adata,
                                       n_top_genes=2000,
                                       flavor="seurat")

    ################################ run method ###############################
    print('\nRunning BANKSY')

    adata.obsm['coord_xy'] = np.array(adata.obs[['x', 'y']])

    if tech.name == '10xVisium':
        n_neighs = 18
    elif tech.name == 'ST':
        n_neighs = 12
    else:
        n_neighs = 15

    banksy_dict = initialize_banksy(adata,
                                    coord_keys=('x', 'y', 'coord_xy'),
                                    num_neighbours=n_neighs,
                                    plt_edge_hist=False,
                                    plt_nbr_weights=False,
                                    plt_theta=False)
    
    banksy_dict, banksy_matrix = generate_banksy_matrix(adata=adata,
                                                        banksy_dict=banksy_dict,
                                                        lambda_list=[0.8],
                                                        max_m=1)

    pca_umap(banksy_dict)
    results_df, _ = run_mclust_partition(banksy_dict=banksy_dict,
                                         num_labels=n_domains)

    ############################### save output ###############################

    labels = results_df.squeeze()['labels']
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
