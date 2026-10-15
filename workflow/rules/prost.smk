rule prost_clustering:
  input:
    os.path.join('workflow','scripts','prost_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','PROST','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','prost_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','PROST','{sample}','prost.log')
  benchmark:
    os.path.join('results','{dataset}','PROST','{sample}','benchmark.tsv')
  conda:
    '../envs/prost.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
