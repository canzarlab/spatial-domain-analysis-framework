rule tacco_clustering:
  input:
    os.path.join('workflow','scripts','tacco_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','TACCO','{sample}','clustering.tsv'),
  params:
    tacco_script = os.path.join('workflow','scripts','tacco_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','TACCO','{sample}','tacco.log')
  benchmark:
    os.path.join('results','{dataset}','TACCO','{sample}','benchmark.tsv')
  conda:
    '../envs/tacco.yml'
  shell:
    'python {params.tacco_script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
    '|| true ; touch {output}'
