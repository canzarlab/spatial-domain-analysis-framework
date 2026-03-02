import glob
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', type=str, required=True)
args = parser.parse_args()

datasets_df = pd.read_csv(os.path.join('config', 'datasets.csv'))
samples_df = pd.read_csv(os.path.join('config', 'samples.csv'))

# find samples associated with the dataset
samplefiles = glob.glob(os.path.join('resources', args.dataset, '*_labels.csv'))
samples = [fl.split('/')[-1].rsplit('_', 1)[0] for fl in samplefiles if 'str_labels' not in fl]

# find samples that are newly added
newsamples = list(set(samples)-set(samples_df.loc[samples_df['dataset_id']==args.dataset, 'sample_id']))
print(f'adding {len(newsamples)} new samples:')
print(newsamples)

# find n_domains associated to each new sample
newndomains = []
for newsample in newsamples:
    labels = pd.read_csv(os.path.join('resources', args.dataset, f'{newsample}_labels.csv'), header=None)
    try:
        uniquelabels = labels.nunique().item()
    except ValueError:
        print(f'there seems to be more than one column in the labels for sample {newsample}')
        raise
    newndomains.append(uniquelabels)

# write to files
if args.dataset not in datasets_df['dataset'].tolist():
    datasets_df.loc[len(datasets_df)] = [args.dataset] + [True]*(len(datasets_df.columns)-1) # untested

for newsample, ndomains in zip(newsamples, newndomains):
    samples_df.loc[len(samples_df)] = [args.dataset, newsample, ndomains] + [True]*(len(samples_df.columns)-3)

datasets_df.to_csv(os.path.join('config', 'datasets.csv'), index=False)
samples_df.to_csv(os.path.join('config', 'samples.csv'), index=False)
