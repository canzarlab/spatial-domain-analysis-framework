rule mln2svg_clustering:
  input:
    os.path.join('workflow','scripts','mln2svg_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','MLN2SVG','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','mln2svg_run.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','MLN2SVG','{sample}','mln2svg.log')
  benchmark:
    os.path.join('results','{dataset}','MLN2SVG','{sample}','benchmark.tsv')
  conda:
    '../envs/mln2svg.yml'
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
