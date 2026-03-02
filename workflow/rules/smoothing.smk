rule smoothing:
  input: 
    os.path.join('results','{dataset}','{method}','{sample}','clustering.tsv'),
    os.path.join('workflow', 'scripts', 'smoothing.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','{method}-smooth','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow', 'scripts', 'smoothing.py'),
  log:
    os.path.join('results', '{dataset}', '{method}-smooth', '{sample}', 'smoothing.log')
  benchmark:
    os.path.join('results','{dataset}','{method}-smooth','{sample}','benchmark.tsv')
  conda:
    '../envs/smoothing.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -m {wildcards.method}'
    '&> {log}'
#    '|| true ; touch {output}'
