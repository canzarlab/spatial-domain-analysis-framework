rule meringue_clustering:
  input:
    os.path.join('workflow','scripts','meringue_run.R'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','MERINGUE','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','meringue_run.R'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','MERINGUE','{sample}','meringue.log')
  benchmark:
    os.path.join('results','{dataset}','MERINGUE','{sample}','benchmark.tsv')
  conda:
    '../envs/meringue.yml'
  shell:
    'Rscript {params.script} {wildcards.dataset} {wildcards.sample} {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
