import os
os.environ["OMP_NUM_THREADS"] = "1"
import pandas as pd
import argparse
import sys; sys.path.insert(0, '/gcm-lfs1/alice/workspace/projSDB')
from workflow.scripts.aggregate import load_coords, load_clustering, unify_labels
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.cluster import adjusted_rand_score, \
                                    fowlkes_mallows_score, \
                                    normalized_mutual_info_score, \
                                    adjusted_mutual_info_score, \
                                    silhouette_score, \
                                    davies_bouldin_score
from sklearn.metrics import balanced_accuracy_score, \
                            accuracy_score, \
                            f1_score, \
                            jaccard_score, \
                            homogeneity_completeness_v_measure
import subprocess
from tqdm import tqdm
from itertools import combinations


def get_ACC(labels_true, labels_pred):
    map_Y = dict(list(zip(pd.unique(labels_true).astype('str').tolist(),
                          list(range(pd.unique(labels_true).size)))))
    map_p = dict(list(zip(pd.unique(labels_pred).astype('str').tolist(),
                          list(range(pd.unique(labels_pred).size)))))
    Y = labels_true.astype('str').map(map_Y).reset_index(drop=True)
    Y_pred = labels_pred.astype('str').map(map_p).reset_index(drop=True)
    assert Y_pred.size == Y.size
    D = max(pd.unique(Y_pred).shape[0], pd.unique(Y).shape[0])+1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(Y_pred.size):
        w[Y_pred[i], Y[i]] += 1
    ind = linear_sum_assignment(w.max() - w)
    acc = sum([w[i, j] for i, j in zip(ind[0], ind[1])])*1.0/Y_pred.size
    return acc


def get_CHAOS(points, labels):
    dist_per_cluster = []
    for k in labels.unique():
        points_cluster = points.loc[labels[labels == k].index]
        if len(points_cluster) < 2:
            continue
        else:
            distances = kneighbors_graph(points_cluster,
                                         n_neighbors=1,
                                         mode='distance')
            dist_per_cluster.append(np.sum(distances))
    chaos = np.sum(dist_per_cluster)/len(labels)
    knns = kneighbors_graph(points, n_neighbors=1, mode='distance') # my addition
    avg_dist = np.median(knns.data)                                 # my addition
    return chaos/avg_dist-1                                         # my addition


def get_PAS(points, labels, n_neighs):
    neighbours = kneighbors_graph(points,
                                  n_neighbors=n_neighs,
                                  mode='connectivity')
    abnormal_list = []
    for ii in range(len(labels)):
        spot_neighbours = neighbours[ii].nonzero()[1]
        neighbour_domains = labels.loc[spot_neighbours]
        spot_domains = labels.loc[ii]
        same_domains = (neighbour_domains == spot_domains).value_counts()
        try:
            if same_domains.loc[False] > n_neighs/2:
                abnormal_list.append(ii)
        except KeyError:
            continue  # means that all neighbours have same label
    return len(abnormal_list)/len(labels)


def get_confusion(labels_df):
    valuecounts = labels_df.apply(lambda col: labels_df.value_counts(['gt', col.name])).fillna(0)
    normvaluecounts = valuecounts.div(valuecounts['gt'].groupby(axis='index', level='gt').sum(), axis='index', level='gt')
    normvaluecounts.drop(columns=['gt'], inplace=True)
    nvc = normvaluecounts.squeeze()
    confusion_dict = {}
    for lbl1, lbl2 in combinations(labels_df['gt'].unique(), 2):
        confusion_dict[(lbl1, lbl2)] = (nvc.get((lbl1, lbl2), default=0) + nvc.get((lbl2, lbl1), default=0)) / (nvc.get((lbl1, lbl1), default=0) + nvc.get((lbl2, lbl2), default=0))
    confusion_df = pd.DataFrame(confusion_dict, index=['confusion']).transpose()
    return confusion_df


def create_and_preprocess_adata_for_umap(X, obs, save_umap=False, savepath_umap=""):
    import scanpy as sc
    adata = sc.AnnData(X, obs)
    sc.pp.filter_genes(adata, min_counts=10)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.neighbors(adata)
    sc.pp.pca(adata)
    sc.tl.umap(adata)
    if save_umap:
        umap_df = pd.DataFrame(adata.obsm['X_umap'])
        if savepath_umap == "" or ".csv" not in savepath_umap:
            raise ValueError("savepath_umap must be specified (including file name, as csv) for UMAP to be saved")
        else:
            umap_df.to_csv(savepath_umap, index=False, header=False)
    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-ds', '--dataset', type=str)
    parser.add_argument('-mt', '--method', type=str)
    parser.add_argument('-sp', '--sample', type=str)
    args = parser.parse_args()

    coords_df = load_coords(ds=args.dataset, sp=args.sample)
    labels_df_not_unified = load_clustering(ds=args.dataset, mt=args.method, sp=args.sample)
    labels_df = unify_labels(labels_df=labels_df_not_unified)

    umap_file = os.path.join('resources', args.dataset, f'{args.sample}_umap.csv')
    try:
        umap_df = pd.read_csv(umap_file, header=None)
    except FileNotFoundError:
        print(f'File {umap_file} does not exist. ')
        print('Creating file...')
        import scanpy as sc
        counts_file = os.path.join('resources', args.dataset, f'{args.sample}_UMImatrix.csv')
        n_lines = int(subprocess.run(['wc', '-l', counts_file],
                                     capture_output=True,
                                     text=True
                                     ).stdout.split(' ')[0])
        counts_df = pd.concat([chunk for chunk in
                               tqdm(pd.read_csv(counts_file, header=None, chunksize=100),
                                    desc='loading count matrix',
                                    total=round(n_lines/100))])
        _ = create_and_preprocess_adata_for_umap(X=counts_df, 
                                                 obs=coords_df, 
                                                 save_umap=True, 
                                                 savepath_umap=umap_file)
        print('File created.')
        umap_df = pd.read_csv(umap_file, header=None)

    if len(labels_df) != len(coords_df):
        print(f'clustering and coords file not the same length \
                ({len(labels_df)} vs {len(coords_df)})')
        failure_file = os.path.join('results',
                                    args.dataset,
                                    args.method,
                                    args.sample,
                                    'cells_deleted.failure')
        with open(failure_file, 'w+') as f:
            f.write('method ran, but deleted cells - fix!')
        raise SystemExit

    labels_noNA_df = labels_df.dropna()
    coords_noNA_df = coords_df.loc[labels_noNA_df.index].reset_index(drop=True)
    umap_noNA_df = umap_df.loc[labels_noNA_df.index].reset_index(drop=True)
    labels_noNA_df = labels_noNA_df.reset_index(drop=True)

    print('calculating supervised metrics')
    ARI = adjusted_rand_score(labels_true=labels_noNA_df['gt'],
                              labels_pred=labels_noNA_df[args.method])
    FM = fowlkes_mallows_score(labels_true=labels_noNA_df['gt'],
                               labels_pred=labels_noNA_df[args.method])
    NMI = normalized_mutual_info_score(labels_true=labels_noNA_df['gt'],
                                       labels_pred=labels_noNA_df[args.method])
    AMI = adjusted_mutual_info_score(labels_true=labels_noNA_df['gt'],
                                     labels_pred=labels_noNA_df[args.method])
    ACC = get_ACC(labels_true=labels_noNA_df['gt'],
                  labels_pred=labels_noNA_df[args.method])
    AC = accuracy_score(y_true=labels_noNA_df['gt'],
                        y_pred=labels_noNA_df[args.method])
    BAC = balanced_accuracy_score(y_true=labels_noNA_df['gt'],
                                  y_pred=labels_noNA_df[args.method])
    FS = f1_score(y_true=labels_noNA_df['gt'],
                  y_pred=labels_noNA_df[args.method],
                  average='weighted')
    JS = jaccard_score(y_true=labels_noNA_df['gt'],
                       y_pred=labels_noNA_df[args.method],
                       average='weighted')
    HOM, COM, VM = homogeneity_completeness_v_measure(labels_true=labels_noNA_df['gt'],
                                                      labels_pred=labels_noNA_df[args.method])

    print('calculating unsupervised metrics')
    SI = silhouette_score(X=coords_noNA_df, labels=labels_noNA_df[args.method])
    DBS = davies_bouldin_score(X=coords_noNA_df, labels=labels_noNA_df[args.method])
    CHAOS = get_CHAOS(points=coords_noNA_df, labels=labels_noNA_df[args.method])
    PAS = get_PAS(points=coords_noNA_df, labels=labels_noNA_df[args.method], n_neighs=10)
    RAS = get_PAS(points=coords_noNA_df, labels=labels_noNA_df[args.method], n_neighs=40)
    SI_X = silhouette_score(X=umap_noNA_df, labels=labels_noNA_df[args.method])
    DBS_X = davies_bouldin_score(X=umap_noNA_df, labels=labels_noNA_df[args.method])
    CHAOS_X = get_CHAOS(points=umap_noNA_df, labels=labels_noNA_df[args.method])
    PAS_X = get_PAS(points=umap_noNA_df, labels=labels_noNA_df[args.method], n_neighs=10)
    RAS_X = get_PAS(points=umap_noNA_df, labels=labels_noNA_df[args.method], n_neighs=40)

    print(f'saving scores to file results/{args.dataset}/{args.method}/{args.sample}/scores.csv')
    scores_df = pd.DataFrame({'dataset': args.dataset,
                              'sample': args.sample,
                              'method': args.method,
                              'ARI': ARI,
                              'FM': FM,
                              'NMI': NMI,
                              'AMI': AMI,
                              'ACC': ACC,
                              'AC': AC,
                              'BAC': BAC,
                              'FS': FS,
                              'JS': JS,
                              'HOM': HOM,
                              'COM': COM,
                              'VM': VM,
                              'SI': SI,
                              'DBS': DBS,
                              'CHAOS': CHAOS,
                              'PAS': PAS,
                              'RAS': RAS,
                              'SI_X': SI_X,
                              'DBS_X': DBS_X,
                              'CHAOS_X': CHAOS_X,
                              'PAS_X': PAS_X,
                              'RAS_X': RAS_X,
                              'num_NA': len(labels_df)-len(labels_noNA_df)},
                             index=['scores'])

    scores_df.to_csv(os.path.join('results',
                                  args.dataset,
                                  args.method,
                                  args.sample,
                                  'scores.csv'),
                     index=False)

    # ensure ground truth scores exist
    if not os.path.exists(os.path.join('results',
                                       args.dataset,
                                       'gt',
                                       args.sample,
                                       'scores.csv')):
        os.makedirs(os.path.join('results', args.dataset, 'gt', args.sample), exist_ok=True)
        print('calculating ground truth scores (unsupervised) and saving them')
        SI = silhouette_score(X=coords_df, labels=labels_df['gt'])
        DBS = davies_bouldin_score(X=coords_df, labels=labels_df['gt'])
        CHAOS = get_CHAOS(points=coords_df, labels=labels_df['gt'])
        PAS = get_PAS(points=coords_df, labels=labels_df['gt'], n_neighs=10)
        RAS = get_PAS(points=coords_df, labels=labels_df['gt'], n_neighs=40)
        SI_X = silhouette_score(X=umap_df, labels=labels_df['gt'])
        DBS_X = davies_bouldin_score(X=umap_df, labels=labels_df['gt'])
        CHAOS_X = get_CHAOS(points=umap_df, labels=labels_df['gt'])
        PAS_X = get_PAS(points=umap_df, labels=labels_df['gt'], n_neighs=10)
        RAS_X = get_PAS(points=umap_df, labels=labels_df['gt'], n_neighs=40)

        scores_df = pd.DataFrame({'dataset': args.dataset,
                                  'sample': args.sample,
                                  'method': 'gt',
                                  'SI': SI,
                                  'DBS': DBS,
                                  'CHAOS': CHAOS,
                                  'PAS': PAS,
                                  'RAS': RAS,
                                  'SI_X': SI_X,
                                  'DBS_X': DBS_X,
                                  'CHAOS_X': CHAOS_X,
                                  'PAS_X': PAS_X,
                                  'RAS_X': RAS_X},
                                 index=['scores'])
        scores_df.to_csv(os.path.join('results',
                                      args.dataset,
                                      'gt',
                                      args.sample,
                                      'scores.csv'),
                         index=False)

    print('done!')
