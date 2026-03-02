rule mnmst_clustering:
  input: 
    os.path.join('workflow','scripts','mnmst_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','MNMST','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','mnmst_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','MNMST','{sample}','mnmst.log')
  benchmark:
    os.path.join('results','{dataset}','MNMST','{sample}','benchmark.tsv')
  conda:
    '../envs/mnmst.yml'
  priority: -1
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
