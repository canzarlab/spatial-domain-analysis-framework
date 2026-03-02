# paper does not actually test on transcriptome data, only imaging

import os
from read_input import read_input

from utag import utag
import scanpy as sc


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'UTAG', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts.to_numpy())
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    lo, hi = 0., 5.
    while hi - lo > 1e-3:
        mi = (lo+hi)/2
        adata = utag(
            adata,
            slide_key=None,
            clustering_method='leiden',
            resolutions=[mi]
        )

        if adata.obs[f'UTAG Label_leiden_{mi}'].nunique() == n_domains:
            break
        if adata.obs[f'UTAG Label_leiden_{mi}'].nunique() < n_domains:
            lo = mi
        else:
            hi = mi

    ############################### save output ###############################
    labels = adata.obs[[f'UTAG Label_leiden_{mi}']]
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
