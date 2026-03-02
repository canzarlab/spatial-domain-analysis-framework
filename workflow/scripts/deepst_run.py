import faulthandler
faulthandler.enable()
import os
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import torch
torch.manual_seed(88)
torch.cuda.manual_seed(88)
import scanpy as sc

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'DeepST/deepst'))
import DeepST


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'DeepST', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    deepen = DeepST.run(save_path=output_folder,
                        task="Identify_Domain",
                        pre_epochs=800,
                        epochs=1000,
                        use_gpu=True)
    spa_typ = 'BallTree'  # LinearRegress not applicable without morphology
    adata = deepen._get_augment(adata, spatial_type=spa_typ, use_morphological=False)
    graph_dict = deepen._get_graph(adata.obsm["spatial"], distType="BallTree")

    # Enhanced data preprocessing
    data = deepen._data_process(adata, pca_n_comps=min(200, adata.n_vars-1, adata.n_obs-1))

    # Training models
    deepst_embed = deepen._fit(data=data, graph_dict=graph_dict,)
    # DeepST outputs
    adata.obsm["DeepST_embed"] = deepst_embed

    adata = deepen._get_cluster_data(adata, n_domains=n_domains, priori=True)

    ############################### save output ###############################

    labels = adata.obs['DeepST_refine_domain']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
