import os
from read_input import read_input
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

import GraphPCA as sg
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances as pair


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'GraphPCA', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    PCA_components = 50

    sc.pp.filter_genes(adata,min_cells=20)
    sc.experimental.pp.normalize_pearson_residuals(adata)
    sc.pp.scale(adata)
    try:
        sc.pp.highly_variable_genes(adata,
#                                    flavor="seurat_v3",
                                    n_top_genes=3000)
    except ValueError:
        sc.pp.highly_variable_genes(adata,
                                    flavor="seurat_v3",
                                    n_top_genes=3000,
                                    span=0.5)

    adata = adata[:, adata.var['highly_variable']]

    shape = 'hexagon' if tech.name == '10xVisium' else 'square' if tech.name == 'ST' else 'generic'
    n_neighs = 6 if shape=='hexagon' else 4 if shape=='square' else 7
    
    location=coords.to_numpy()

    Z,W = sg.Run_GPCA(adata, location=location, n_components = 50, method = "knn", platform = "Visium", _lambda =0.5,
                      save_reconstruction=True)
    adata.obsm["GraphPCA"] = Z

    estimator = KMeans(n_clusters=n_domains)
    res = estimator.fit(Z[:,:])
    lable_pred=res.labels_
    adata.obs["GPCA_pred"]= lable_pred
    adata.obs["GPCA_pred"] = adata.obs["GPCA_pred"].astype('category')
    refined_pred=sg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["GPCA_pred"].tolist(), dis=pair(location), shape=shape, neighbor_num=n_neighs)
    adata.obs["GPCA_pred"]= refined_pred
    adata.obs["GPCA_pred"] = adata.obs["GPCA_pred"].astype('category')
    
    ############################### save output ###############################

    labels = adata.obs[['GPCA_pred']]
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
