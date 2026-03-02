rule graphpca_clustering:
  input:
    os.path.join('workflow','scripts','graphpca_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','GraphPCA','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','graphpca_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','GraphPCA','{sample}','graphpca.log')
  benchmark:
    os.path.join('results','{dataset}','GraphPCA','{sample}','benchmark.tsv')
  conda:
    '../envs/graphpca.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
