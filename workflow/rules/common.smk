rule get_n_bins:
    input:
        os.path.join('resources', '{dataset}', '{sample}_coords.csv'),
        script = os.path.join('workflow', 'scripts', 'bin_numbers.py'),
    output:
        os.path.join('resources', '{dataset}', '{sample}_nbins.csv'),
    log:
        os.path.join('resources', '{dataset}', '{sample}_nbins.log'),
    conda:
        '../envs/basics.yml',
    shell:
        'python {input.script} -ds {wildcards.dataset} -sp {wildcards.sample}'
        '&> {log}'


def get_n_domains(wildcards):
    n_domains = sp_df[(sp_df.dataset_id == wildcards.dataset) & (sp_df.sample_id == wildcards.sample)].n_cluster.values[0]
    return n_domains


rule get_scores:
    input:
        os.path.join('results', '{dataset}', '{method}', '{sample}', 'clustering.tsv'),
        script = os.path.join('workflow', 'scripts', 'get_scores.py'),
    output: 
        os.path.join('results', '{dataset}', '{method}', '{sample}', 'scores.csv'),
    log:
        os.path.join('results', '{dataset}', '{method}', '{sample}', 'scores.log'),
    conda:
        '../envs/basics.yml',
    priority: 2
    shell:
        'python {input.script} -ds {wildcards.dataset} -mt {wildcards.method} -sp {wildcards.sample}'
        '&> {log}'


rule collect_scores:
    input: 
        all_scores_files = [os.path.join('results',
                                         row.dataset_id,
                                         row.method_id,
                                         row.sample_id,
                                         'scores.csv')
                            for row in sp_df.itertuples()],
    output: 
        file = os.path.join('results', 'all_scores.csv'),
    run:
        all_scores = []
        for scores_file in input.all_scores_files:
            all_scores.append(pd.read_csv(scores_file))
        pd.concat(all_scores).to_csv(output.file, index=False)
