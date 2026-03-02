import pandas as pd
import os
import glob
from tqdm import tqdm


def number_of_optimal_bins(coords, nbins_col, nbins_row):
    indexcol = pd.cut(coords['col'], nbins_col, labels=range(nbins_col))
    indexrow = pd.cut(coords['row'], nbins_row, labels=range(nbins_row))
    index_coords = pd.DataFrame({'index_col': indexcol,
                                 'index_row': indexrow}).astype('int')
    n_duplicates = pd.DataFrame({'duplicates': index_coords.value_counts()})
    n_neighbours = pd.DataFrame(index_coords.groupby(['index_col', 'index_row'])
                                            .size()
                                            .unstack()
                                            .pipe(lambda df: df.where(df.isna(), 1)).fillna(0)
                                            .pipe(lambda df: df*(df.shift(-1, axis=0).fillna(0)
                                                                 + df.shift(1, axis=0).fillna(0)
                                                                 + df.shift(-1, axis=1).fillna(0)
                                                                 + df.shift(1, axis=1).fillna(0)))
                                            .stack()
                                            .loc[index_coords.drop_duplicates().set_index(['index_col', 'index_row']).index],
                                index=pd.MultiIndex.from_frame(index_coords.drop_duplicates()),
                                columns=['neighbours'])
    vcs = n_duplicates.join(n_neighbours).value_counts()
    try:
        return (True, vcs.loc[(1, 4.0)] + vcs.loc[(1, 3.0)])
    except KeyError:
        return (False, vcs.reset_index()[['duplicates', 'neighbours']].mean())


def get_optimal_values(coords, colrange, rowrange):
    print(f'getting optimal values for colrange {colrange} and rowrange {rowrange}')
    maxes = []
    values = []
    for ii in colrange:
        for jj in rowrange:
            passed, ret = number_of_optimal_bins(coords, ii, jj)
            if passed:
                maxes.append([ii, jj, ret])
            else:
                values.append([ii, jj, *ret.values.tolist()])
    if len(maxes) > 0:
        df = pd.DataFrame(maxes, columns=['col', 'row', 'n_optimal']).sort_values('n_optimal')
        if colrange.step == 1 and rowrange.step == 1:
            retval = df.iloc[-1][['col', 'row']].values.tolist()
            return retval
        df_best = df.iloc[-1-min(5, len(df)-1):-1][['col', 'row']]
        colmin, rowmin = df_best.min()
        colmax, rowmax = df_best.max()
        colmax = max(colmin+1, colmax)
        rowmax = max(rowmin+1, rowmax)
        colstep = max(round((colmax - colmin)/20), 1)
        rowstep = max(round((rowmax - rowmin)/20), 1)
        retval = get_optimal_values(coords,
                                    range(colmin, colmax, colstep),
                                    range(rowmin, rowmax, rowstep))
        return retval
    else:
        raise UserWarning
        df = pd.DataFrame(values, columns=['col', 'row', 'n_dupes', 'n_neighs'])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', '-ds', type=str, required=True)
    parser.add_argument('--sample', '-sp', type=str, required=True)
    args = parser.parse_args()

    dataset = args.dataset
    sample = args.sample
    fname = os.path.join('resources', dataset, f'{sample}_coords.csv')
    coords = pd.read_csv(fname, header=None, names=['col', 'row'])

    ll = len(coords)
    if ll > 100000:
        colmin, colmax, colstep = round(ll/10000), round(ll/100), round(ll/500)
        rowmin, rowmax, rowstep = round(ll/10000), round(ll/100), round(ll/500)
    elif ll > 10000:
        colmin, colmax, colstep = round(ll/1000), round(ll/10), round(ll/200)
        rowmin, rowmax, rowstep = round(ll/1000), round(ll/10), round(ll/200)
    else:
        colmin, colmax, colstep = round(ll/200), round(ll/10), round(ll/100)
        rowmin, rowmax, rowstep = round(ll/200), round(ll/10), round(ll/100)

    optcol, optrow = get_optimal_values(coords,
                                        range(colmin, colmax, colstep),
                                        range(rowmin, rowmax, rowstep))

    nbins_df = pd.DataFrame({'sample': [sample], 'row': [optrow], 'col': [optcol]})
    nbins_df.to_csv(os.path.join('resources', dataset, f'{sample}_nbins.csv'), index=False)
