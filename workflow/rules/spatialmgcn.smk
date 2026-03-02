rule spatialmgcn_clustering:
  input: 
    os.path.join('workflow','scripts','spatialmgcn_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','SpatialMGCN','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','spatialmgcn_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','SpatialMGCN','{sample}','spatialmgcn.log')
  benchmark:
    os.path.join('results','{dataset}','SpatialMGCN','{sample}','benchmark.tsv')
  conda:
    '../envs/spatialmgcn.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
