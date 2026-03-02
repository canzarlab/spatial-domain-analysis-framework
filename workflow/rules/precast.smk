rule precast_clustering:
  input:
    os.path.join('workflow','scripts','precast_run.R'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','PRECAST','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','precast_run.R'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','PRECAST','{sample}','precast.log')
  benchmark:
    os.path.join('results','{dataset}','PRECAST','{sample}','benchmark.tsv')
  conda:
    '../envs/precast.yml'
  shell:
    'Rscript {params.script} {wildcards.dataset} {wildcards.sample} {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
