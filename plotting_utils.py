import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from typing import Optional, Union
import tqdm
from workflow.scripts.aggregate import load_scores


def load_and_annotate_scores(scores_to_include='all', include_gt=False):
    all_scores_df = load_scores()
    all_scores_df['technology'] = [ds.split('_')[0] for ds in all_scores_df['dataset']]
    all_scores_df['tech_type'] = ['imaging' if tech in ['MERFISH','osmFISH'] else 'high-res' if tech in ['Slide-seq','STARmap'] else 'low-res' if tech in ['10xVisium','ST'] else 'synth' for tech in all_scores_df['technology']]
    all_scores_df['exp'] = [sp.split('-')[0] if 'exp' in sp else None for sp in all_scores_df['sample']]

    if not include_gt:
        all_scores_df = all_scores_df[all_scores_df['method']!='gt']

    real_datasets = ['osmFISH_Codeluppi', 'MERFISH_Moffitt', 'MERFISH_Zhang', 'STARmap_Wang', 'Slide-seq_Langlieb', '10xVisium_Fu', '10xVisium_Maynard', 'ST_Stahl']
    synth_datasets = ['synth-close', 'synth-circles', 'synth-moons', 'synth-varied', 'synth-blobs', 'synth-trials', 'synth-2blobs', 'synth-3stripes']

    if scores_to_include == 'real':
        all_scores_df = all_scores_df[all_scores_df['dataset'].isin(real_datasets)]
        return all_scores_df
    elif scores_to_include == 'synth':
        all_scores_df = all_scores_df[all_scores_df['dataset'].isin(synth_datasets)]
        return all_scores_df
    elif scores_to_include == 'all':
        all_scores_df['type'] = all_scores_df['dataset'].apply(lambda f: 'real' if f in real_datasets else 'synth' if f in synth_datasets else 'other')
        return all_scores_df
    else: 
        raise ValueError(f'{scores_to_include} not known: choose from "real", "synth" or "all".')
    

def plot_heatmap_fig2(plotdf):
    fig, axes = plt.subplots(ncols=2, figsize=(8.3, 2.8), gridspec_kw={'width_ratios': (28, 1), 'wspace': 0.05}, layout='constrained')
    vmin = -0.4
    sns.heatmap(plotdf.transpose(), 
                cmap='coolwarm', cbar_kws={'label': r'ARI - max(ARI$_{base}$)', 'orientation': 'vertical'}, cbar_ax=axes[1],
                center=0, 
                vmin=vmin,
                annot=plotdf.fillna('').mask(plotdf<vmin).fillna('<').mask(plotdf>=vmin).fillna('').transpose(), fmt='',
                square=True, ax=axes[0])
    axes[0].set_xlabel('')
    axes[0].set_ylabel('')
    return fig
    

def convert_pvalue_to_asterisks(pvalue):
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"


def plot_merfish_visium_fig2(scores_df):
    plotdf = scores_df[scores_df['dataset'].isin(['MERFISH_Zhang','10xVisium_Maynard'])]
    order = plotdf.pivot_table(values='ARI', index='method', columns='technology', aggfunc='median')
    order['diff'] = order['10xVisium']-order['MERFISH']
    order = order.sort_values('diff').index
    fig, axes = plt.subplots(nrows=2, figsize=(10, 4), gridspec_kw={'height_ratios': (1, 4), 'hspace': 0.1})
    pas_df = plotdf.pivot_table(values='PAS', index='method', columns='technology', aggfunc='mean')
    pas_df['mean'] = pas_df.mean(axis='columns')
    
    sns.barplot(plotdf, x='method', y='PAS', estimator='mean', hue='technology', order=order, 
                palette={'MERFISH': "#c53e84", '10xVisium': "#9bb964"}, 
                width=0.6, capsize=.1, err_kws={'linewidth': 1}, legend=False, ax=axes[0])
    padding = 0.8
    axes[0].set_xlim(0-padding, 29+padding)
    axes[0].set_xticks([])
    axes[0].set_xlabel('')

    sns.boxplot(plotdf, x='method', y='ARI', hue='technology', 
                palette={'MERFISH': "#c53e84", '10xVisium': "#9bb964"}, gap=.25, fliersize=2, order=order, ax=axes[1])

    axes[1].get_legend().remove()
    handles, labels = axes[1].get_legend_handles_labels()
    axes[1].legend(handles=handles, labels=['MERFISH\N{EN DASH}Zhang', 'Visium\N{EN DASH}Maynard'], loc='upper right', bbox_to_anchor=(.6, -0.5))

    plt.tick_params('x', labelrotation=90)
    axes[1].set_ylim(-0.1, 0.98)
    axes[1].set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    axes[1].set_xlabel('')

    def get_better_on_dict_merfish_visium():
        plotdf = scores_df[scores_df['dataset'].isin(['MERFISH_Zhang','10xVisium_Maynard'])]
        better_on_dict = {}
        for mt in plotdf['method'].unique():
            aris_visium = plotdf[(plotdf['technology']=='10xVisium')&(plotdf['method']==mt)]['ARI']
            aris_merfish = plotdf[(plotdf['technology']=='MERFISH')&(plotdf['method']==mt)]['ARI']
            pvalue = mannwhitneyu(aris_visium, aris_merfish).pvalue
            mark = convert_pvalue_to_asterisks(pvalue)
            if '*' in mark and aris_visium.median() > aris_merfish.median():
                better_on_dict[mt] = (mark, 'significantly better on visium')
            elif '*' in mark and aris_visium.median() <= aris_merfish.median():
                better_on_dict[mt] = (mark, 'significantly better on merfish')
            else:
                better_on_dict[mt] = (mark, 'not significantly different')
        return better_on_dict
        
    merfishvisiumpalette = {'significantly better on merfish': "#c53e84", 'not significantly different': '#8f8f8f', 'significantly better on visium': "#9bb964"}
    better_on_dict = get_better_on_dict_merfish_visium()
    for ii, mt in enumerate(order):
        try:
            x1, x2 = ii-0.3, ii+0.3
            y, h, col = 0.85, 0.02, 'k'
            axes[1].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
            axes[1].text((x1+x2)/2, y+h, better_on_dict[mt][0], fontsize=9, ha='center', va='bottom', color=col)
            axes[1].axvline(ii, color=merfishvisiumpalette[better_on_dict[mt][1]], linewidth=17.5, alpha=.2)
            axes[1].set_xlim([-.75, len(order)-0.25])
        except ValueError: 
            continue


def plot_stability_fig2(scores_df):
    plotdf = scores_df[(scores_df['dataset'].str.contains('_perturbed'))&(scores_df['sample'].str.contains('init'))&(~scores_df['sample'].str.contains('maynard'))&(~scores_df['sample'].str.contains('MERFISH'))&(scores_df['method']!='gt')]
    plotdf = plotdf[~plotdf['method'].isin(['scanpy','Seurat','scanpy-smooth','Seurat-smooth'])]
    plotdf[['orig_sample', 'pert', 'seed']] = [sp.rsplit('_', 2) for sp in plotdf['sample']]
    plotdf = plotdf.sort_values(['dataset','method','orig_sample','seed'])
    method_goal_median = plotdf.groupby(['dataset','method'])['ARI'].median()
    method_sample_median = plotdf.groupby(['dataset','method','orig_sample'])['ARI'].median()
    method_sample_mediandifference = method_sample_median-method_goal_median
    plotdf['ARI (sample bias removed)'] = [row.ARI-method_sample_mediandifference.loc[row.dataset, row.method, row.orig_sample] for row in plotdf.itertuples(index=False)]
    plotdf = plotdf.set_index('dataset')
    fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(8, 3), gridspec_kw={'height_ratios': (2, 3), 'hspace': 0}, layout='constrained')
    stds_df = plotdf.groupby(['dataset','method'])['ARI (sample bias removed)'].std().reset_index()
    order = stds_df.groupby('method')['ARI (sample bias removed)'].mean().sort_values().index
    hue_order = ['MERFISH_Zhang_perturbed','10xVisium_Maynard_perturbed']
    palette = ["#c53e84", "#9bb964"]
    sns.barplot(stds_df, x='method', y='ARI (sample bias removed)', hue='dataset', order=order, hue_order=hue_order, palette=palette, alpha=.8, ax=axes[0])
    axes[0].get_legend().remove()
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, ['MERFISH\N{EN DASH}Zhang', 'Visium\N{EN DASH}Maynard'], bbox_to_anchor=(0, .9), loc='center left', frameon=False)
    sns.boxplot(plotdf, hue='dataset', x='method', y='ARI (sample bias removed)', order=order, hue_order=hue_order, 
                fliersize=1, saturation=.6, palette=palette, linecolor='dimgray', linewidth=.8, fill=False, legend=False, ax=axes[1])
    for ax in [axes[0]]:
        ax.set_xlabel('')
        ax.set_ylabel('s.d.')
        ax.tick_params('x', bottom=False, labelbottom=False)
        ax.spines['top'].set_visible(False)
    for ax in [axes[1]]:
        ax.tick_params('x', labelrotation=90)
        ax.set_xlabel('')
        ax.set_ylabel('ARI')
        ax.set_ylim(-0.1, 0.85)


def load_synth_scores():
    synth_scores_df = load_and_annotate_scores('synth')
    synthdf = synth_scores_df[~synth_scores_df['sample'].str.contains('exp2-labels')][['sample','method','ARI','PAS','RAS','SI_X','technology','exp']]
    synthdf['shape'] = ['stripes' if 'close'in tech else tech.split('-')[1] for tech in synthdf['technology']]
    synthdf['sample'] = [sp.rsplit('-', 1)[0].split('-', 1)[1] for sp in synthdf['sample']]
    synthdf = synthdf[(~synthdf['sample'].str.contains('-sd1'))&(synthdf['shape']!='trials')]
    synthdf['div'] = [row.sample[3:] if row.exp=='exp1' else np.nan for row in synthdf.itertuples(index=False)]
    synthdf.loc[:, ('exp2label', 'exp2prop')] = [(row.sample.split('-')[0], np.float64(row.sample.split('-')[1][4:])) if row.exp=='exp2' else (np.nan, np.nan) for row in synthdf.itertuples(index=False)]
    synthdf.loc[:, ('exp3type', 'exp3label', 'exp3prop')] = [(row.sample.split('-')[0], row.sample.split('-')[1], np.float64(row.sample.split('-')[2][4:])) if row.exp=='exp3' else (np.nan, np.nan, np.nan) for row in synthdf.itertuples(index=False)]
    synthdf['exp4prop'] = [np.float64(row.sample.split('-')[1][4:]) if row.exp=='exp4' else np.nan for row in synthdf.itertuples(index=False)]
    synthdf['exp5prop'] = [np.float64(row.sample[4:]) if row.exp=='exp5' else np.nan for row in synthdf.itertuples(index=False)]
    synthdf['exp6prop'] = [np.float64(row.sample[9:]) if row.exp=='exp6' else np.nan for row in synthdf.itertuples(index=False)]
    synthdf['exp7type'] = [row.sample[4:9] if row.exp=='exp7' else np.nan for row in synthdf.itertuples(index=False)]
    synthdf['exp7prop'] = [np.float64(row.sample[9:]) if row.exp=='exp7' else np.nan for row in synthdf.itertuples(index=False)]
    synthdf['exp8length'] = [np.float64(row.sample[10:]) if row.exp=='exp8' else np.nan for row in synthdf.itertuples(index=False)]
    synthdf['exp9perm'] = [row.sample if row.exp=='exp9' else np.nan for row in synthdf.itertuples(index=False)]

    synthdf['exp2label'] = ['label1or2' if lbl in ['label1','label2'] else lbl for lbl in synthdf['exp2label']]
    synthdf = synthdf[(~synthdf['exp2prop'].isin([0.25, 0.75, 0.95, 0.99]))]
    synthdf['exp3type'] = ['add2or1' if lbl in ['add1','add2'] else lbl for lbl in synthdf['exp3type']]
    synthdf['exp3label'] = ['label1or2' if lbl in ['label1','label2'] else lbl for lbl in synthdf['exp3label']]
    synthdf = synthdf[(~synthdf['exp3prop'].isin([0.25, 0.75]))]

    synthdf.drop(columns=['technology', 'sample'], inplace=True)
    synthdf.sort_values(['exp','method','shape','div'], inplace=True)
    return synthdf


def plot_thickness_fscores_fig4(fscores_df):
    exp1celltypecolours = {1: '#3355ff', 2: '#9944ff',3: '#ffaa00'}
    palette = [exp1celltypecolours[2], exp1celltypecolours[3]]
    def middle_domain_conditions(row):
        return (row.expname=='1, 2 similar' and row.celltype=='ct2') or (row.expname=='1, 3 similar' and row.celltype=='ct3')
    def right_domain_conditions(row):
        return (row.expname=='1, 2 similar' and row.celltype=='ct3') or (row.expname=='1, 3 similar' and row.celltype=='ct2')
    fscores_df['thickness'] = fscores_df.index.to_frame().apply(lambda row: row['div']-50 if middle_domain_conditions(row) else 100-row['div'] if right_domain_conditions(row) else 50, axis='columns')
    fscores_df = fscores_df.loc[:, :, :, ['ct2', 'ct3'], :]
    mtgroups = [['BASS','ADEPT','TACCO','SpiceMix','SpaceFlow','SCAN-IT'],['DeepST','SpatialPCA','STAGATE','PAST','MNMST'],['SpaDo','CellCharter','UTAG','Vesalius'],['BayesSpace','SC-MEB','GraphST','SEDR'],['GraphPCA','SpaGCN','SpatialMGCN','BANKSY','MERINGUE'],['Seurat-smooth','Seurat','PRECAST','scanpy','CCST','scanpy-smooth']]
    fig, axes = plt.subplots(ncols=3, nrows=2, sharex=True, sharey=True, figsize=(7, 4), layout='constrained', gridspec_kw={'hspace': 0.1})
    for ii, (mtgroup, ax) in tqdm.tqdm(enumerate(zip(mtgroups, axes.reshape(-1))), total=len(mtgroups)):
        sns.lineplot(fscores_df.loc[:, :, :, :, mtgroup], x='thickness', y='fscore', hue='celltype', palette=palette, estimator='mean', ax=ax)
        ax.set_xticks([5, 25, 45])
        ax.set_title(f"group {ii+1}")
        ax.get_legend().remove()
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles=handles, labels=["similar" if lbl=="ct2" else "dissimilar" for lbl in labels], bbox_to_anchor=(1, 0.5), loc='center left', fontsize='medium', title='cell type')


def plot_shapedependence_fig4(synthdf):
    plotdf = synthdf[synthdf['exp']=='exp9']
    plotdf = plotdf[plotdf['exp9perm'].map(len)==3]
    order = plotdf.groupby(['method','shape'])['ARI'].mean().groupby('method').median().sort_values(ascending=False).index

    plotdf['ARI-mean(ARI)'] = (plotdf.set_index('method')['ARI']-plotdf.groupby('method')['ARI'].mean()).tolist()

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 5), gridspec_kw={'height_ratios': (2, 4), 'width_ratios': (50, 1)}, layout='constrained')
    sns.barplot(plotdf, x='method', y='ARI', estimator='mean', order=order, color='gray', ax=axes[0, 0])
    axes[0, 0].tick_params('x', labelrotation=90)
    axes[0, 1].remove()
    axes[0, 0].set_ylabel('mean ARI')
    plotdf = plotdf.pivot_table(values='ARI-mean(ARI)', index='shape', columns='method')[order]
    vmax = 0.11
    sns.heatmap(plotdf, cmap='PRGn', center=0, vmax=vmax, vmin=-vmax, annot=plotdf.round(2).mask((plotdf<vmax)&(plotdf>-vmax)).fillna(''), fmt='', cbar_ax=axes[1, 1], cbar_kws={'label': 'ARI difference to mean'}, ax=axes[1, 0])
    axes[1, 0].tick_params('y', rotation=0)
    axes[1, 0].set_yticklabels(['circular', 'laminar', 'concentric', 'complex'])
    for ax in axes.reshape(-1):
        ax.set_xlabel('')


def plot_spots(
        coords_df: pd.DataFrame,
        labels_df: Union[pd.DataFrame, pd.Series],
        colours: Optional[Union[np.ndarray, tuple, dict]] = None,
        markersize: Optional[np.float64] = 0, 
        rasterized: Optional[bool] = True,
        ax: Optional[plt.Axes] = None,
        save_path: Optional[str] = None) -> plt.Axes:
    '''
    plots spots with given labels.

    Args:
        coords_df (pd.DataFrame): "x" and "y" coordinates to plot
        labels_df (Union[pd.DataFrame, pd.Series]): labels for each spot,
                                    header will be plot title
        colours (Optional[np.ndarray]): list of colours for each domain
        markersize (Optional[np.float64]): size of markers
        rasterized (Optional[bool]): whether to rasterize plot
        ax (Optional[plt.Axes]): ax to plot in
        save_path (Optional[str]): name and path under which to save plot

    Returns:
        plt.Axes: ax with plot
    '''
    if not ax:
        _, ax = plt.subplots(1, 1)

    if not set(coords_df.columns) == {'x', 'y'}:
        raise AssertionError('expect columns "x" and "y" in coords_df')

    if not coords_df.shape[0] == labels_df.shape[0]:
        raise AssertionError('coords_df and labels_df not the same size ' +
                             f'({coords_df.shape[0]} vs {labels_df.shape[0]})')

    if not all(coords_df.index == labels_df.index):
        raise AssertionError('index of coords_df not the same as labels_df, ' +
                             'unable to assign labels to spots (same shape)')

    if isinstance(labels_df, pd.DataFrame):
        if len(labels_df.columns) > 1:
            raise AssertionError('labels_df has too many columns: ' +
                                 f'{len(labels_df.columns)} (should have 1)')
        labels = labels_df[labels_df.columns[0]]
    if isinstance(labels_df, pd.Series):
        labels = labels_df

    if isinstance(colours, np.ndarray) and len(colours) < len(labels.unique()):
        raise AssertionError(f'n_colours ({len(colours)}) smaller than ' +
                             f'n_labels ({len(labels.unique())})')

    if isinstance(colours, dict) and len(set(labels.unique()) - set(colours.keys())) > 0:
        raise AssertionError('n_colours dict is not complete.')

    if markersize == 0:
        markersize = max(6-labels.shape[0]/1000, 0.5)

    coords_df_scaled = (coords_df - coords_df.min()) / \
                       (coords_df.max() - coords_df.min())

    labels_noNA = [label for label in labels.unique() if not pd.isna(label)]
    try:
        labels_noNA = [int(label) for label in labels_noNA]
    except ValueError:
        print('labels are not (castable to) integers')

    for ii, label in enumerate(labels_noNA):
        label_index = labels[labels == label].index
        coords_subset_df = coords_df_scaled.loc[label_index]
        if isinstance(colours, np.ndarray) or isinstance(colours, tuple):
            ax.plot(coords_subset_df['x'], coords_subset_df['y'],
                    '.', markersize=markersize,
                    rasterized=rasterized, # set to False for better quality
                    color=colours[ii])
        elif isinstance(colours, dict):
            ax.plot(coords_subset_df['x'], coords_subset_df['y'],
                    '.', markersize=markersize,
                    rasterized=rasterized,
                    color=colours[label])
        else:
            try:
                import colorcet as cc
                colors = cc.glasbey_category10
            except ModuleNotFoundError:
                print("module 'colorcet' not available, use default colours")
                colors = plt.cm.tab10.colors
            try: # if this doesn't work, the colors will not match
                color = colors[label]
            except TypeError:
                print(f'since label {label} for {labels.name} is not int, colours will be random')
                color = colors[ii]
            ax.plot(coords_subset_df['x'], coords_subset_df['y'],
                    '.', markersize=markersize,
                    rasterized=rasterized,
                    color=color)

    # ax.set_aspect('equal', 'box')
    ax.set_title(labels.name, wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])

    if save_path:
        plt.savefig(save_path, dpi=600)


def get_plotsize(coords_df: pd.DataFrame) -> tuple[int, int]:
    if not set(coords_df.columns) == {'x', 'y'}:
        raise AssertionError('expect columns "x" and "y" in coords_df')

    ranges = coords_df.max() - coords_df.min()
    return (np.sqrt(6*ranges['y']/ranges['x']),
            np.sqrt(6*ranges['x']/ranges['y']))


def plot_spots_multilabels(
        coords_df: pd.DataFrame,
        labels_df: pd.DataFrame,
        title: Optional[str] = None,
        unify: Optional[bool] = True,
        colours: Optional[Union[np.ndarray, tuple, dict]] = None,
        max_ncols: Optional[int] = 6,
        save_path: Optional[str] = None,
        ratio_mod: Optional[np.float64] = 1) -> plt.Figure:
    '''
    plots spots with multiple sets of labels.

    Args:
        coords_df (pd.DataFrame): "x" and "y" coordinates to plot
        labels_df (pd.DataFrame): columns containing sets of labels
        title (Optional[str]): title
        unify (Optional[bool]): whether to unify the cluster labels between methods
        colours (Optional[np.ndarray]): list of colours for each domain
        max_ncols (Optional[int]): maximal number of columns, default: 6
        save_path (Optional[str]): name under which to save plot
        ratio_mod (Optional[np.float64]): ratio modifier for plot axes

    Returns:
        plt.Figure: figure with one plot for each set of labels
    '''
    if len(labels_df.columns) == 1:
        print('labels_df has only one column, use plot_spots')

    nlabels = len(labels_df.columns)
    scale_x, scale_y = get_plotsize(coords_df)
    scale_x = scale_x/ratio_mod
    scale_y = scale_y*ratio_mod
    ncols = min(nlabels, max_ncols)
    nrows = nlabels//max_ncols+1

    fig = plt.figure(figsize=(scale_x * ncols,
                              scale_y * (nlabels//max_ncols+1)))

    if unify:
        from workflow.scripts.aggregate import unify_labels
        labels_df = unify_labels(labels_df)

    for ii, labels in enumerate(labels_df.columns):
        ax = plt.subplot(nrows, ncols, ii+1)
        plot_spots(coords_df=coords_df,
                   labels_df=labels_df[labels],
                   colours=colours, 
                   ax=ax)

    if title:
        fig.suptitle(title)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, bbox_inches='tight', dpi=600)


def plot_spots_multisamples(df_dict: dict[str,
                                          tuple[pd.DataFrame,
                                                pd.DataFrame]],
                            title: Optional[str] = None,
                            max_ncols: Optional[int] = 6,
                            save_path: Optional[str] = None) -> plt.Figure:
    '''
    plots labels of different samples in one plot.

    Args:
        df_dict (dict): dictionary with subplot titles as keys and
                        (coords_df, labels_df) tuples as values
        title (Optional[str]): title
        max_ncols (Optional[int]): maximal number of columns, default: 6
        save_path (Optional[str]): name under which to save plot

    Returns:
        plt.Figure: figure with one plot for each sample
    '''
    nsamples = len(df_dict.keys())
    ncols = min(nsamples, max_ncols)
    nrows = nsamples//max_ncols+1

    fig = plt.figure(figsize=(3*ncols, 3*nrows))
    for ii, (name, (coords, labels)) in enumerate(df_dict.items()):
        ax = plt.subplot(nrows, ncols, ii+1)
        labels.name = name
        plot_spots(coords_df=coords,
                   labels_df=labels,
                   ax=ax)

    if title:
        fig.suptitle(title)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, bbox_inches='tight', dpi=600)
