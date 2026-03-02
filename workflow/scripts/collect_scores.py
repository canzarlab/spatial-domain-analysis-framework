import os
import pandas as pd
from aggregate import load_designmatrix
import tqdm

sp_df = load_designmatrix()

all_scores = []
no_scores = []
no_gt_score = []
for row in tqdm.tqdm(sp_df.itertuples(index=False), total=len(sp_df)):
    try:
        scores = pd.read_csv(os.path.join('results',
                                          row.dataset_id,
                                          row.method_id,
                                          row.sample_id,
                                          'scores.csv'))
        all_scores.append(scores)
    except FileNotFoundError:
        no_scores.append(row)
        continue

ds_sp_df = sp_df[['dataset_id','sample_id']].drop_duplicates()
for row in tqdm.tqdm(ds_sp_df.itertuples(index=False), total=len(ds_sp_df)):
    try:
        scores = pd.read_csv(os.path.join('results',
                                          row.dataset_id,
                                          'gt',
                                          row.sample_id,
                                          'scores.csv'))
        all_scores.append(scores)
    except:
        no_gt_score.append(row)
        continue

all_scores_df = pd.concat(all_scores).drop_duplicates()
all_scores_df.to_csv(os.path.join('results', 'scores.csv'), index=False)

print(f'scores files found for {100*(1-len(no_scores)/len(sp_df)):.2f}% of targets and {100*(1-len(no_gt_score)/len(ds_sp_df)):.2f}% of ground truths')