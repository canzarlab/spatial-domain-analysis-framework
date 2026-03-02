rule banksy_clustering:
  input: 
    os.path.join('workflow','scripts','banksy_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','BANKSY','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','banksy_run.py'),
    n_domains = get_n_domains,
  log:    
    os.path.join('results','{dataset}','BANKSY','{sample}','banksy.log')
  benchmark:
    os.path.join('results','{dataset}','BANKSY','{sample}','benchmark.tsv')
  conda:
    '../envs/banksy.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
