import os
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

import tacco as tc

from read_input import read_input


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'TACCO', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    adata = sc.AnnData(X=counts, obs=coords)

    ################################ run method ###############################
    print('\nRunning method')

    adata = tc.preprocessing.filter(adata, min_counts_per_gene=10, min_counts_per_cell=10,
                            assume_valid_counts=True)  # otherwise fails for non-integer
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    print('calculating regions')
    lo, hi = -20., 20.
    while hi - lo > 1e-3:
        print(f'lo: {lo}, hi: {hi}')
        mi = (lo + hi)/2
        tc.tl.find_regions(adata, resolution=mi)
        n_pred = len(adata.obs['region'].unique())
        if n_pred == n_domains:
            lo = mi
            break
        elif n_pred < n_domains:
            lo = mi
        else:
            hi = mi
    print(f'{sample} leiden res: {lo:.2f}')

    ############################### save output ###############################
    labels = adata.obs['region']
    try:
        assert len(labels) == len(coords)
    except AssertionError as e:
        print(e)
        print('inserting back missing cells.')
        missing_indices = set([f'{ii}' for ii in range(len(coords))])-set(labels.index)
        for idx in missing_indices:
            labels.loc[idx] = ''
        assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
