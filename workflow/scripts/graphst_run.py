import os
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import torch
torch.manual_seed(88)
torch.cuda.manual_seed(88)
import scanpy as sc

from GraphST import GraphST
from GraphST.utils import clustering


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'GraphST', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')
    GraphST.preprocess(adata)
    GraphST.construct_interaction(adata)

    # datatype is left at default "10x" for all except Slide-seq data
    datatype = 'Slide' if 'Slide' in tech.name else '10x'
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model = GraphST.GraphST(adata,
                            device=device,
                            datatype=datatype)
    adata = model.train()

    clustering(adata=adata,
               n_clusters=n_domains)

    ############################### save output ###############################

    labels = adata.obs['domain']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
