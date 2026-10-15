rule sctm_clustering:
  input:
    os.path.join('workflow','scripts','sctm_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','scTM','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','sctm_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','scTM','{sample}','sctm.log')
  benchmark:
    os.path.join('results','{dataset}','scTM','{sample}','benchmark.tsv')
  conda:
    '../envs/sctm.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'