import os
import argparse
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

from scipy.spatial.distance import pdist, squareform


def refine(pred, adj_df, nbs):
    refined = []
    for ii, self_pred in enumerate(pred):
        adj_tmp = adj_df.loc[ii, :].sort_values()
        nbs_df = adj_tmp[0:nbs+1]
        nbs_pred = pred.loc[nbs_df.index]
        label_count = nbs_pred.value_counts()

        if self_pred not in label_count.index or label_count.loc[self_pred] <= nbs/2:
            refined.append(label_count.idxmax())
        else:
            refined.append(self_pred)
    return pd.DataFrame(refined)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', '-d', type=str, required=True)
    parser.add_argument('--sample', '-s', type=str, required=True)
    parser.add_argument('--method', '-m', type=str, required=True)
    parser.add_argument('--technology', '-t', type=str)
    args = parser.parse_args()
    dataset = args.dataset
    sample = args.sample
    method = args.method
    tech = args.technology if args.technology else dataset.split('_')[0]
    tech = pd.read_csv(os.path.join('config', 'technology_meta.csv'),
    				   index_col=0).loc[tech]
    assert isinstance(int(tech.n_neighbours), int)

    ################################ read data ################################
    print('Reading data')

    input_folder = os.path.join('resources', dataset)
    orig_folder = os.path.join('results', dataset, method, sample)
    output_folder = os.path.join('results', dataset, f'{method}-smooth', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  
    coords = pd.read_csv(os.path.join(input_folder, f'{sample}_coords.csv'), 
     					 header=None, 
     					 names=['x', 'y'])
    orig_labels = pd.read_csv(os.path.join(orig_folder, 'clustering.tsv'),
                              header=0, 
                              names=['label']
                              ).squeeze()
    assert orig_labels.shape[0] == coords.shape[0]

    ################################ run method ###############################
    print('\nRunning method')

    adjacency = pd.DataFrame(squareform(pdist(coords)))

    n_nbs = int(tech.n_neighbours)

    labels = refine(pred=orig_labels,
                    adj_df=adjacency,
                    nbs=n_nbs)

    ############################### save output ###############################

    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
