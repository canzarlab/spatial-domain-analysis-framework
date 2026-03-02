rule scmeb_clustering:
  input:
    os.path.join('workflow','scripts','scmeb_run.R'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
    os.path.join('resources','{dataset}','{sample}_nbins.csv'),
  output:
    os.path.join('results','{dataset}','SC-MEB','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','scmeb_run.R'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','SC-MEB','{sample}','scmeb.log')
  benchmark:
    os.path.join('results','{dataset}','SC-MEB','{sample}','benchmark.tsv')
  conda:
    '../envs/scmeb.yml'
  shell:
    'Rscript {params.script} {wildcards.dataset} {wildcards.sample} {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
