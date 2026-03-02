rule utag_clustering:
  input: 
    os.path.join('workflow','scripts','utag_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','UTAG','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','utag_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','UTAG','{sample}','utag.log')
  benchmark:
    os.path.join('results','{dataset}','UTAG','{sample}','benchmark.tsv')
  conda:
    '../envs/utag.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
