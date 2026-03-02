import os
from read_input import read_input
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import squidpy as sq
import scipy as sp
import numpy as np



if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'SpiceMix', sample, 'files')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  

    adata = sc.AnnData(X=counts, obs=coords)
    
    ############################# preprocess data #############################
    print('\nPreprocessing\n')
    
    coords.to_csv(os.path.join(output_folder, f'meta_{sample}.csv'), 
                  sep=',', 
                  index=None)

    sc.pp.filter_genes(adata, min_cells=round(adata.n_obs/10))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    counts = pd.DataFrame(adata.X)
    counts.to_csv(os.path.join(output_folder, f'expression_{sample}.txt'), 
                  sep='\t', 
                  header=None, 
                  index=None)

    genes = pd.DataFrame(counts.index)
    genes.to_csv(os.path.join(output_folder, f'genes_{sample}.txt'),
                 header=None,
                 index=None)

    adata.obsm['spatial'] = adata.obs.to_numpy()
    sq.gr.spatial_neighbors(adata, delaunay=True)

    indices = sp.sparse.find(adata.obsp['spatial_connectivities'])
    edges = np.concatenate((indices[0], indices[1])).reshape((-1,2), order='F')
    edges = pd.DataFrame(edges)
    edges.to_csv(os.path.join(output_folder, f'neighborhood_{sample}.txt'), 
                 sep='\t', 
                 header=None, 
                 index=None)