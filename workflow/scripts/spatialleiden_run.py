import os
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import numpy as np
np.random.seed(88)
import scanpy as sc
import squidpy as sq
import spatialleiden as sl

RANDOM_STATE = 88


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'SpatialLeiden', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    print('Preprocessing')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, random_state=RANDOM_STATE)
    sc.pp.neighbors(adata, random_state=RANDOM_STATE)

    print('Building spatial graph')
    sq.gr.spatial_neighbors(adata,
                            coord_type='generic',
                            n_neighs=tech.n_neighbours)
    adata.obsp['spatial_connectivities'] = sl.distance2connectivity(
        adata.obsp['spatial_distances'])

    print('Searching resolutions for target cluster count')
    spatial_kwargs = {'layer_ratio': 1.8,
                      'random_state': RANDOM_STATE,
                      'directed': (False, True)}
    latent_resolution, spatial_resolution = sl.search_resolution(
        adata,
        n_domains,
        latent_kwargs={'random_state': RANDOM_STATE},
        spatial_kwargs=spatial_kwargs)

    print('Clustering')
    sl.spatialleiden(adata,
                     resolution=(latent_resolution, spatial_resolution),
                     layer_ratio=1.8,
                     directed=(False, True),
                     random_state=RANDOM_STATE)

    ############################### save output ###############################

    labels = adata.obs['spatialleiden']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
