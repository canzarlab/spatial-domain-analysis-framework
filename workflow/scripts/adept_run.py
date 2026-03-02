import os
import argparse
from tqdm import tqdm
import subprocess
import pandas as pd
import numpy as np
np.random.seed(88)
import scanpy as sc
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'ADEPT/'))

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import os
from imputation.impute import impute_
import GAAE as GAAE
import argparse
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
from GAAE.utils import impute, DE_num_calc

import sklearn


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', '-d', type=str, required=True)
    parser.add_argument('--sample', '-s', type=str, required=True)
    parser.add_argument('--n_domains', '-n', type=int, required=True)
    parser.add_argument('--technology', '-t', type=str)
    parser.add_argument('--data_dir', type=str, default="./", help="root dir for input data")
    parser.add_argument('--gt_dir', type=str, default="./", help="root dir for data ground truth")
    parser.add_argument('--input_data', type=str, default="151673", help="input data section id")
    parser.add_argument('--impute_cluster_num', type=str, default="7", help="diff cluster numbers for imputation")
    parser.add_argument('--cluster_num', type=int, default=7, help="input data cluster number")
    parser.add_argument('--radius', type=int, default=150, help="input data radius")
    parser.add_argument('--no_de', type=int, default=0, help='switch on/off DEG selection module')
    parser.add_argument("--use_mean", type=int, default=0, help="use mean value in de list or not")
    parser.add_argument("--impute_runs", type=int, default=1, help="time of runs for imputation")
    parser.add_argument("--runs", type=int, default=2, help="total runs for the data")
    parser.add_argument('--use_hvgs', type=int, default=3000, help="select highly variable genes before training")
    parser.add_argument('--use_preprocessing', type=int, default=1, help='use preprocessed input or raw input')
    parser.add_argument('--save_fig', type=int, default=0, help='saving output visualization')
    parser.add_argument('--de_nzr_min', type=float, default=0.299, help='suggested min nzr threshold after de selection')
    parser.add_argument('--de_nzr_max', type=float, default=0.399, help='suggested max nzr threshold after de selection')
    parser.add_argument('--use_gpu_id', type=str, default='0', help='use which GPU, only applies when you have multiple gpu')

    args = parser.parse_args()
    dataset = args.dataset
    sample = args.sample
    n_domains = args.n_domains

    args.impute_cluster_num = [str(n_domains)]
    args.filter_nzr = 0.15
    args.filter_num = 5
    root_d = args.data_dir

    
    tech = args.technology if args.technology else dataset.split('_')[0]
    tech = pd.read_csv(os.path.join('config', 'technology_meta.csv'),
                       index_col=0).loc[tech]
    assert tech.seq_or_img in ['seq', 'img']

    ################################ read data ################################
    print('Reading data')

    input_folder = os.path.join('resources', dataset)
    output_folder = os.path.join('results', dataset, 'ADEPT', sample)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    file = os.path.join(input_folder, f'{sample}_UMImatrix.csv')
    n_lines = int(subprocess.run(['wc', '-l', file],
                                 capture_output=True,
                                 text=True
                                 ).stdout.split(' ')[0])
    counts = pd.concat([chunk for chunk
                        in tqdm(pd.read_csv(file, header=None, chunksize=100),
                                desc='count matrix progress',
                                total=round(n_lines/100))])
    coords = pd.read_csv(os.path.join(input_folder, f'{sample}_coords.csv'),
                         header=None,
                         names=['x', 'y'])
    assert counts.shape[0] == coords.shape[0]
    adata = sc.AnnData(X=counts)
    coords = coords.to_numpy()
    adata.obsm['spatial'] = coords
    
    lo, hi = 0, 400
    knn = 6
    BT = sklearn.neighbors.BallTree(coords)
    while hi - lo > 1e-3:
        mi = (lo + hi)/2
        cnt = BT.query_radius(coords, mi, count_only=True)
        if np.mean(cnt) > knn: hi = mi
        else: lo = mi

    args.radius = lo
    ################################ run method ###############################
    print('\nRunning method')
    adata.obs['original_clusters'] = np.ones(len(adata.X))
    adata.obs['Ground Truth'] = np.ones(len(adata.X))

    gene_min_count = 5
    sc.pp.filter_genes(adata, min_counts=gene_min_count)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    de_top_k_list = [] #DE_num_calc(args, adata)

    de_list_epoch = []
    if de_top_k_list != []:
        print("performing DEGs selection")
        adata_list = []
        for de_ in de_top_k_list:
            for cluster_n in args.impute_cluster_num: # not working
                print("cluster_n = ", cluster_n)
                GAAE.get_kNN(adata, rad_cutoff=args.radius)

                ari_ini, ari_final, de_list, adata_out = GAAE.train_ADEPT_use_DE(adata, n_epochs=1000,
                                                                                 num_cluster=int(cluster_n),
                                                                                 dif_k=de_, device_id=args.use_gpu_id)
                de_list_epoch.append(de_list)
                adata_list.append(adata_out)
                        
        g_union = [] if len(de_list_epoch) == 0 else set.union(*de_list_epoch)
        imputed_ad = impute(args, adata_list, g_union, de_top_k_list)
    else:
        print("skip performing DEGs selection")
        imputed_ad = adata

    """result of imputed data"""
    GAAE.get_kNN(imputed_ad, rad_cutoff=args.radius)
    ari_ini, ARI, de_list, adata_out = GAAE.train_ADEPT_use_DE(imputed_ad, n_epochs=1000, num_cluster=n_domains, device_id=args.use_gpu_id)


    ############################### save output ###############################

    labels = adata_out.obs['mclust_impute']
    assert len(labels) == len(coords)

    labels.to_csv(os.path.join(output_folder, 'clustering.tsv'),
                  sep='\t',
                  header=['region'],
                  index=False)
    print('\n Saved domain labels to clustering.tsv')
