import os
from read_input import read_input
import torch

import squidpy as sq
import cellcharter as cc
import scanpy as sc
import scvi
import numpy as np

if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'CellCharter', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts.to_numpy())
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    scvi.settings.seed = 12345
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
    )
    model = scvi.model.SCVI(adata)
    model.train(early_stopping=True, enable_progress_bar=True)
    adata.obsm['X_scVI'] = model.get_latent_representation(adata)\
                                .astype(np.float32)

    sq.gr.spatial_neighbors(adata,
                            coord_type='generic',
                            delaunay=True,
                            spatial_key='spatial')
    cc.gr.remove_long_links(adata)
    cc.gr.aggregate_neighbors(adata,
                              n_layers=3,
                              use_rep='X_scVI')
    if torch.cuda.is_available():
        trainer_params = dict(accelerator='gpu', devices=1)
    else:
        trainer_params = None
    autok = cc.tl.ClusterAutoK(
        n_clusters=(n_domains, n_domains),
        model_params=dict(
            random_state=12345,
            # If running on GPU
            trainer_params=trainer_params
        )
    )
    autok.fit(adata, use_rep='X_cellcharter')
    adata.obs['cluster_cellcharter'] = autok.predict(adata,
                                                     use_rep='X_cellcharter')

    ############################### save output ###############################
    labels = adata.obs[['cluster_cellcharter']]
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
