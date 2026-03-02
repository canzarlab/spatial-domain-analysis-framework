import os
from read_input import read_input
os.environ['OMP_NUM_THREADS'] = '2'
os.environ['MKL_NUM_THREADS'] = '2'
os.environ['NUMEXPR_NUM_THREADS'] = '2'
import scanpy as sc
import faulthandler
faulthandler.enable()
import warnings
warnings.filterwarnings('ignore')

from SpaceFlow import SpaceFlow


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'SpaceFlow', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    sc.pp.filter_genes(adata, min_counts=1)  # necessary for preprocessing
    sf = SpaceFlow.SpaceFlow(adata=adata)
    sf.preprocessing_data(n_top_genes=3000)

    sf.train(embedding_save_filepath=os.path.join(output_folder, 'embedding.tsv'))

    lo, hi = 0., 5.
    while hi - lo > 1e-4:
        mi = (lo + hi)/2
        sf.segmentation(domain_label_save_filepath=os.path.join(output_folder, 'domains.tsv'),
                        resolution=mi)
        n_pred = len(sf.domains.unique())
        if n_pred == n_domains:
            lo = mi
            break
        elif n_pred < n_domains:
            lo = mi
        else:
            hi = mi
    print(f'segmentation res: {lo:.2f}')

    ############################### save output ###############################

    labels = sf.domains
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
