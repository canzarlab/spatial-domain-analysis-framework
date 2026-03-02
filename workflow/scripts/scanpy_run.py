import os
from read_input import read_input
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'scanpy', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts, obs=coords)

    ################################ run method ###############################
    print('\nRunning method')

    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.neighbors(adata)

    lo, hi = -20., 20.
    while hi - lo > 1e-3:
        print(f'lo: {lo}, hi: {hi}')
        mi = (lo + hi)/2
        sc.tl.leiden(adata, resolution=mi)
        n_pred = len(adata.obs['leiden'].unique())
        if n_pred == n_domains:
            lo = mi
            break
        elif n_pred < n_domains:
            lo = mi
        else:
            hi = mi
    print(f'{sample} leiden res: {lo:.2f}')

    ############################### save output ###############################

    labels = adata.obs['leiden']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
