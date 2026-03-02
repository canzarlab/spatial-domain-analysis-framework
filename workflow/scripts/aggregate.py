import os
import pandas as pd
import numpy as np
from typing import Optional
from collections import defaultdict
import warnings


def load_designmatrix(load_all=False) -> pd.DataFrame:
    '''
    load dataframe of all sp/mt combinations to run

    Args:
        load_all (bool): whether to load all combinations,
                         irrespective of "include" status

    Returns:
        pd.DataFrame: ds/sp/mt combinations with n_cluster and platform
    '''
    datasets_df = pd.read_csv('config/datasets.csv', sep=',')
    datasets_included = datasets_df.loc[datasets_df.include]['dataset']
    methods_df = pd.read_csv('config/methods.csv', sep=',')
    methods_included = methods_df.loc[methods_df.include]['method']
    samples_df = pd.read_csv('config/samples.csv', sep=',')
    samples_df = samples_df.fillna('')

    # create design matrix
    sp_df = pd.melt(
        samples_df,
        # id_vars=['dataset_id', 'sample_id', 'n_cluster', 'platform'],
        id_vars=['dataset_id', 'sample_id', 'n_cluster'],
        var_name='method_id', value_name='include'
        )
    # remove unused datasets, methods and sp/mt combos
    if not load_all:
        sp_df = sp_df[sp_df.include].drop(columns='include')

    if not load_all:
        sp_df = sp_df.loc[sp_df['dataset_id'].isin(datasets_included)]
        sp_df = sp_df.loc[sp_df['method_id'].isin(methods_included)]

    return sp_df


def load_scores(
        file: Optional[str] = 'results/scores.csv') -> pd.DataFrame:
    '''
    load scores file to dataframe

    Args:
        file (Optional[str]): path to scores, default results/all_scores.csv

    Returns:
        pd.DataFrame: scores
    '''
    # read in data
    scores = pd.read_csv(file)
    # get columns corresponding to scores
    numcols = scores.drop(columns=['sample', 'method', 'dataset'])
    # get score names
    score_names = numcols.columns.tolist()
    # make score columns numeric and remove nans
    numcols = numcols.apply(pd.to_numeric, errors='coerce').fillna(0)
    scores[score_names] = numcols

    # list of datasets
    # datasets = scores['dataset'].unique()

    # techs = {ds: ds.rsplit('_', 1)[0] for ds in datasets}
    # scores['tech'] = scores['dataset'].map(techs)

    scores = scores.drop_duplicates()

    return scores  # , score_names


def load_coords(
        ds: str,
        sp: str) -> pd.DataFrame:
    '''
    loads coordinate file with x and y columns

    Args:
        ds (str): name of dataset
        sp (str): name of sample

    Returns:
        pd.DataFrame: coordinates
    '''
    filename = os.path.join('resources',
                            ds,
                            f'{sp}_coords.csv')
    try:
        coords_df = pd.read_csv(filename, header=None, names=['x', 'y'])
    except FileNotFoundError:
        print(f'no coords file at {filename}')
        raise
    return coords_df


def load_labels(
        ds: str,
        sp: str) -> pd.DataFrame:
    '''
    loads ground truth labels

    Args:
        ds (str): name of dataset
        sp (str): name of sample

    Returns:
        pd.DataFrame: labels
    '''
    filename = os.path.join('resources', ds, f'{sp}_labels.csv')
    try:
        return pd.read_csv(filename, header=None, names=['gt'])
    except FileNotFoundError:
        warnings.warn(f'no ground truth labels found for {ds}, {sp}')
        return None


def load_clustering(
        ds: str,
        mt: str,
        sp: str,
        include_gt: bool = True) -> pd.DataFrame:
    '''
    loads a given clustering result

    Args:
        ds (str): name of dataset
        mt (str): name of method
        sp (str): name of sample
        include_gt (bool): whether to include "gt" column

    Returns:
        pd.DataFrame: clusterings
    '''
    labels_gt = load_labels(ds=ds, sp=sp) if include_gt else pd.DataFrame()
    filepath = os.path.join('results', ds, mt, sp, 'clustering.tsv')

    if not os.path.exists(filepath):
        raise UserWarning('does not exist')
    elif os.stat(filepath).st_size == 0:
        raise UserWarning('is empty')
    else:
        clust = pd.read_csv(filepath, sep='\t', header=0, names=[mt])
        if clust[mt].isna().all():
            raise UserWarning('contains only NA')
        clust.reset_index(drop=True, inplace=True)

    if include_gt and len(labels_gt) != len(clust):
        warnings.warn('ground truth and clustering are of different lengths')

    return clust.join(labels_gt)


def load_clusterings(
        ds: str,
        sp: str,
        include_gt: bool = True,
        print_warnings: bool = False) -> pd.DataFrame:
    '''
    loads clustering results for given sample

    Args:
        ds (str): name of dataset
        sp (str): name of sample
        include_gt (bool): whether to include "gt" column
        print_warnings (bool): whether to print warnings

    Returns:
        pd.DataFrame: clusterings
    '''
    if include_gt:
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                labels_gt = load_labels(ds=ds, sp=sp)
            except UserWarning:
                include_gt = False

    # add label assignments for each method
    sp_df = load_designmatrix()
    sp_mt_df = sp_df[(sp_df['dataset_id'] == ds) &
                     (sp_df['sample_id'] == sp)]
    labels_list = [labels_gt] if include_gt else []
    failed_dict = defaultdict(list)
    for mt in sp_mt_df['method_id'].unique():
        try:
            clust = load_clustering(ds=ds, mt=mt, sp=sp, include_gt=False)
            labels_list.append(clust)
        except UserWarning as w:
            failed_dict[w].append(mt)

    if print_warnings and len(failed_dict) > 0:
        print('\nproblems with clustering.tsv:')
        [print(f'  {issue} for {", ".join(failed_dict[issue])}')
         for issue in failed_dict.keys()]
        print('these methods will be discarded.')

    set_of_lengths = set([df.shape for df in labels_list])
    if not len(set_of_lengths) == 1:
        if include_gt:
            if print_warnings: 
                print("\nwill proceed with labels of gt shape, but:")
            good_shape = labels_gt.shape
        else:
            if print_warnings:
                print("\nwill proceed with labels of majority shape, but:")
            good_shape = pd.DataFrame([df.shape for df in labels_list])\
                           .value_counts().idxmax()

        problem_list = [(*labels.columns.tolist(), labels.shape[0])
                        for labels in labels_list
                        if labels.shape != good_shape]
        if print_warnings:
            [print(f'  {mt} has {n} labels instead of {good_shape[0]}.')
             for (mt, n) in problem_list]
            print('these methods will be discarded.')
            print('to avoid this, ensure labels for each coordinate.')

        labels_list = [labels for labels in labels_list
                       if labels.shape == good_shape]

    labels_df = pd.concat(labels_list, axis=1)

    return labels_df


def unify_labels(
        labels_df: pd.DataFrame) -> pd.DataFrame:
    '''
    unifies label names within a set of clusterings 

    Unifies to column "gt" if it exists, otherwise to first column

    Args:
        labels_df (pd.DataFrame): input clusterings

    Returns:
        pd.DataFrame: unified clusterings
    '''
    from scipy.optimize import linear_sum_assignment

    def numberfy_labels(labels):
        labels_noNA = [str(label) for label in pd.unique(labels) if not pd.isna(label)]
        map_labels = dict(zip(labels_noNA,
                              list(range(len(labels_noNA)))))
        return labels.astype('str').map(map_labels).reset_index(drop=True)
    labels_numbers = labels_df.apply(numberfy_labels, axis=0)

    if 'gt' in labels_numbers.columns:
        lb_gt = labels_numbers['gt']
    else: 
        lb_gt = labels_numbers.dropna(axis=1).iloc[:,0]

    def get_mapping(lb_pred):
        assert len(lb_gt) == len(lb_pred)
        lb_pred.reset_index(drop=True)
        lb_gt.reset_index(drop=True)
        lb_pred_noNA = lb_pred.dropna().astype('int')
        lb_gt_match = lb_gt.loc[lb_pred_noNA.index].astype('int')
        assert len(lb_pred_noNA) == len(lb_gt_match)
        n_spots = len(lb_gt_match)
        w = np.zeros((max(lb_gt_match)+1, max(lb_pred_noNA)+1),
                     dtype=np.int64)
        for i in range(n_spots):
            w[lb_gt_match.iloc[i], lb_pred_noNA.iloc[i]] += 1
        mapping = linear_sum_assignment(w.max() - w)
        return dict(zip(mapping[1], mapping[0]))
    mappings = labels_numbers.apply(get_mapping, axis=0)
    labels_df_unified = labels_numbers.apply(
                        lambda col: col.map(mappings[col.name]))
    return labels_df_unified


def get_assignments(
        ds: str,
        sp: str,
        aris_in_title: bool = False,
        print_warnings: bool = False):
    '''
    loads coordinates and unified cluster assignments for a given sample

    Args:
        ds (str): name of dataset
        sp (str): name of sample
        aris_in_title (bool): whether to add ARIs in column names
        print_warnings (bool): whether to print warnings

    Returns:
        pd.DataFrame: coordinates
        pd.DataFrame: cluster assignments
    '''
    # load assignments by each method for this specific sample
    coords_df = load_coords(ds=ds, sp=sp)
    labels_df = load_clusterings(ds=ds, sp=sp, print_warnings=print_warnings)

    if not coords_df.shape[0] == labels_df.shape[0]:
        raise ValueError('coords and labels not the same length: ' +
                         f'{coords_df.shape[0]} vs {labels_df.shape[0]}')

    # add title including ARI
    if aris_in_title:
        titles = []
        for mt in labels_df.columns:
            if mt != 'gt':
                try:
                    ari = get_score(ds, sp, mt)['ARI']
                except KeyError:
                    ari = 0
                titles.append(f'{mt}: {ari:.2f}')
            else:
                titles.append(mt)
        # titles = [f'{mt}: {get_score(ds, sp, mt)["ARI"]:.2f}'
        #           if mt != 'gt' else mt
        #           for mt in labels_df.columns]
        labels_df.columns = titles

    # unify assignment labels (i.e. make labels similar to ground truth)
    if 'gt' in labels_df.columns:
        labels_df_unified = unify_labels(labels_df)
    else:
        labels_df_unified = labels_df

    return coords_df, labels_df_unified


def get_adjusted_assignments(
        ds: str,
        sp: str):
    '''
    loads coordinates and labels adjusted for ground truth correspondence

    mostly relevant for synthetic samples!

    Args:
        ds (str)
        sp (str)

    Returns:
        pd.DataFrame: coordinates
        pd.DataFrame: adjusted cluster assignments
    '''
    coords_df, labels_df = get_assignments(ds, sp)
    origlabels = pd.read_csv(f'resources/{ds}/{sp}_labels.csv', header=None, names=['orig'])
    correspondence = labels_df[['gt']].join(origlabels).drop_duplicates().set_index('gt')['orig'].to_dict()
    newlabels = labels_df.apply(lambda col: col.map(correspondence))
    return coords_df, newlabels


def get_adjusted_labels(ds: str,
                        sp: str,
                        mt: str) -> pd.Series:
    '''
    gets labels adjusted for ground truth correspondence (mostly relevant for synthetic samples!)

    Args:
        ds (str)
        sp (str)
        mt (str)

    Returns:
        pd.Series: labels
    '''
    labels_with_gt = load_clustering(ds, mt, sp, include_gt=True)
    labels_unified = unify_labels(labels_with_gt)
    correspondence = labels_with_gt[['gt']].rename(columns={'gt': 'orig'}).join(labels_unified['gt']).drop_duplicates().set_index('gt')['orig'].to_dict()
    adjustedlabels = labels_unified[mt].map(correspondence)
    return adjustedlabels


def get_score(
        dataset: str,
        sample: str,
        method: str,
        scores_df: Optional[pd.DataFrame] = None) -> pd.Series:
    '''
    gets scores for a given ds/sp/mt

    Args:
        dataset (str)
        sample (str)
        method (str)
        scores_df (Optional[pd.DataFrame]): scores file

    Returns:
        pd.Series: scores
    '''
    if not isinstance(scores_df, pd.DataFrame):
        scores_df = load_scores()

    scores_df = scores_df.set_index(['dataset', 'sample', 'method'])

    try:
        scores = scores_df.loc[dataset, sample, method]
    except KeyError:
        raise KeyError(f'no scores for {dataset}/{sample}/{method}')

    return scores
