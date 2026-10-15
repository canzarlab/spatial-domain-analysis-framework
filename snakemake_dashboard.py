import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd
pd.set_option('max_colwidth', None)
import datetime as dt
from tabulate import tabulate
import re

log_dir = '.snakemake/log/'
log_file = os.listdir(log_dir)[-1]
print(f'reading log file {log_file}')

with open(log_dir + log_file, 'r') as log_file:
    log = log_file.readlines()

# filter irrelevant lines
if len(log) == 0:
    print('there is nothing yet!')
    raise SystemExit
else:
    log = [line for line in log if 'Activating conda environment:' not in line
                                and 'Select jobs to execute' not in line
                                and 'Waiting at most 5 seconds' not in line
                                and line[0] != '[']

# break up log into blocks, separated by empty lines
blocks = {}
n_empty = 0
blocks[n_empty] = []
for line in log:
    if line == '\n' or line.split()[0] in ['rule', 'Error', 'MissingOutputException', 'Finished']:
        n_empty += 1
        blocks[n_empty] = [line]
    else:
        blocks[n_empty].append(line.strip())

blocks = [block for block in blocks.values() if len(block) > 0 and block[0]!='\n']
short_blocks = [block for block in blocks if len(block)==1]

# get numbers of total jobs to execute from first block
total_jobs_list = []
try:
    table_start = [nr for nr, line in enumerate(blocks[0]) if '---' in line][0]
except IndexError:
    print(blocks[0][0])
    raise SystemExit
header_row = re.split(r'( ){2,}', blocks[0][table_start-1])[slice(0,8,2)]
for line in blocks[0][table_start+1:]:
    total_jobs_list.append(line.split())
total_jobs_df = pd.DataFrame(total_jobs_list, columns=header_row)
total_jobs_df.set_index('job', inplace=True)

# evaluate rest of the blocks for job status
jobs_dict = {}
for ii, block in enumerate(blocks[1:]):
    try:
        if block[0].split()[0] == 'rule':
            jobid = [ln.split()[-1] for ln in block if 'jobid' in ln][0]
            try:
                jobs_dict[jobid] = dict(zip([block[0].split()[0]] + [ln.split(':')[0].strip() for ln in block[1:]],
                                        [block[0].split()[1].strip(':')] + [ln.split(':')[1] for ln in block[1:]]))
            except IndexError:
                print(block)
            jobs_dict[jobid]['started'] = 1
            jobs_dict[jobid]['failed'] = 0
            jobs_dict[jobid]['missing'] = 0
            jobs_dict[jobid]['done'] = 0
        elif block[0].split()[0] == 'Error':
            jobid = [ln.split()[-1] for ln in block if 'jobid' in ln][0]
            jobs_dict[jobid]['failed'] = 1
        elif block[0].split()[0] == 'MissingOutputException':
            jobid = [ln.split()[-1] for ln in block if 'completed successfully, but' in ln][0]
            jobs_dict[jobid]['missing'] = 1
        elif block[0].split()[0] == 'Finished':
            jobid = block[0].split()[-1].strip('.')
            jobs_dict[jobid]['done'] = 1
        else:
            print(f'new block type: {block[0]}')
    except:
        print('exciting new error in this block: \n\n')
        print(block)
jobs_df = pd.DataFrame(jobs_dict).transpose()
try:
    dropcols = ['jobid', 'resources', 'benchmark', 'input']
    jobs_df = jobs_df.drop(columns=dropcols)
except KeyError as e:
    undroppable = str(e).split('[')[1].split(']')[0]
    jobs_df = jobs_df.drop(columns=[col for col in dropcols if col not in undroppable])

# calculate stuff to print
n_started = len(jobs_df)
n_done = len(jobs_df[jobs_df['done'] == 1])
n_total = int(total_jobs_df.loc['total', 'count'])
start_time = dt.datetime.strptime(log_file.name.split('/')[-1].split('.')[0],
                                  "%Y-%m-%dT%H%M%S")
curr_time = dt.datetime.today()

# print stuff
print(f'{n_started} of {n_total} jobs ({100*n_started/n_total:.2f}%) started and {n_done} ({100*n_done/n_total:.2f}%) done.')
if 'Complete log:' in log[-1]:
    print('\nsnakemake completed!')
else:
    print(f'\nprojected end: {(start_time+((curr_time-start_time)/(n_started/n_total))).strftime("%a %b %d, %H:%M")}\n')

# create dataframe to show
new_df = jobs_df.pivot_table(values=['failed', 'missing', 'done', 'started'], index='rule', aggfunc='sum')
show_df = pd.DataFrame(index=[job for job in total_jobs_df.index if job!='total'])
show_df[new_df.columns] = new_df
show_df['total'] = total_jobs_df['count'].astype(int)
show_df['progress'] = 100*(show_df['done']+show_df['failed']+show_df['missing'])/show_df['total']
show_df['progress'] = show_df['progress'].apply(lambda col: f'{col if col>0 else 0:.1f}%')
show_df = show_df.sort_values('total')
show_df.drop(columns=['missing', 'started'], inplace=True)

# print dataframe
print(
    tabulate(
        show_df.astype('string')
                       .fillna('')
                       # .map(lambda x: x.split('.')[0])
                       # .join(erroreds['most recent failure'])
                       # .fillna('')
        , headers='keys'
        , tablefmt='github'
        , numalign='right'
        , stralign='right'
    )
)

# print errors
if sum(jobs_df['failed']) + sum(jobs_df['missing']) > 0:
    errs_df = pd.melt(jobs_df[(jobs_df['missing'] == 1) | (jobs_df['failed'] == 1)][['rule', 'log', 'failed', 'missing']], id_vars=['rule', 'log'], value_vars=['failed', 'missing'], var_name='outcome')
    errs_df = errs_df[errs_df['value'] == 1].drop(columns=['value'])[['rule', 'outcome', 'log']].sort_values('outcome')
    errs_df.to_csv('errors.csv', index=False)
    errs_df = errs_df.drop_duplicates(['rule', 'outcome'], keep='last').set_index('rule').sort_values('outcome')
    errs_df.drop(columns=['outcome'], inplace=True)
    print('\nproblematic jobs:\n')
    print(tabulate(errs_df, headers='keys', tablefmt='simple'))

# print currently running
if sum(jobs_df['started']) > sum(jobs_df['missing']) + sum(jobs_df['done']) + sum(jobs_df['failed']):
    print('\ncurrently running:')
    running_df = jobs_df[(jobs_df['started'] == 1) & (jobs_df['done'] == 0) & (jobs_df['failed'] == 0) & (jobs_df['missing'] == 0)]
    running_df = running_df[running_df['rule'] != 'get_scores']
    try:
        running_df.drop(columns=['rule', 'output', 'reason', 'wildcards', 'priority', 'started', 'failed', 'missing', 'done'], inplace=True)
    except KeyError:
        running_df.drop(columns=['rule', 'output', 'reason', 'wildcards', 'started', 'failed', 'missing', 'done'], inplace=True)
    running_df.set_index('log', inplace=True)
    print(tabulate(running_df.sort_index(), tablefmt='simple'))
