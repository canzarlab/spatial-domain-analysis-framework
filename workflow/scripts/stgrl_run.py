import os
import sys
from read_input import read_input
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import scanpy as sc
import torch
torch.manual_seed(88)
torch.cuda.manual_seed(88)

# Point rpy2 at the R installation shipped inside the conda env so mclust is
# importable (the upstream example hardcodes an env-specific R_HOME).
os.environ.setdefault('R_HOME', os.path.join(sys.prefix, 'lib', 'R'))

# stGRL uses flat, relative imports (e.g. `from preprocess import ...`), so the
# repository directory has to be on sys.path rather than imported as a package.
sys.path.insert(0, os.path.join('workflow', 'scripts', 'stGRL'))
from stGRL import stGRL
from utils import clustering


if __name__ == '__main__':
    dataset, sample, n_domains, tech, counts, coords = read_input()
    output_folder = os.path.join('results', dataset, 'stGRL', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adata = sc.AnnData(X=counts)
    adata.obsm['spatial'] = coords.to_numpy()

    ################################ run method ###############################
    print('\nRunning method')

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    # Sparse path for imaging/large-spot technologies, dense otherwise.
    datatype = 'Slide' if tech.seq_or_img == 'img' else '10X'

    model = stGRL(adata, device=device, datatype=datatype)
    adata = model.train()

    clustering(adata,
               n_domains,
               radius=tech.n_neighbours,
               method='mclust',
               refinement=True)

    ############################### save output ###############################

    labels = adata.obs['stGRL']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
