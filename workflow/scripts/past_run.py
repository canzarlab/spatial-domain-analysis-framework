import os
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')

import past
import scanpy as sc
import torch

past.setup_seed(88)


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'PAST', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

    print('preprocessing')
    # not using default values for n_tops and gene_method because this is
    # used in all tutorials
    n_tops = min(adata.shape[1], 3000)
    adata = past.preprocess(adata,
                            n_tops=n_tops,
                            gene_method='gearyc')

    print('setting up model')
    d_lat = min(adata.shape[1]-1, 50)
    PAST_model = past.PAST(d_in=adata.shape[1],
                           d_lat=d_lat
                           ).to(device)

    print('training model')
    PAST_model.model_train(adata, device=device)

    print('clustering')
    adata = PAST_model.output(adata)
    adata = past.mclust_R(adata, num_cluster=n_domains)

    ############################### save output ###############################

    labels = adata.obs[['mclust']]
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
