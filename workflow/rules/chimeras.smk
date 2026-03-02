rule chimera1_clustering:
  input: 
    os.path.join('workflow','scripts','chimera_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','chimera-1','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','chimera_run.py'),
    n_domains = get_n_domains,
    outfile = os.path.join('results','{dataset}','chimera-1','{sample}','clustering.tsv'),
    outdir = os.path.join('results','{dataset}','chimera-1'),
  log:
    os.path.join('results','{dataset}','chimera-1','{sample}','chimera1.log'),
  benchmark:
    os.path.join('results','{dataset}','chimera-1','{sample}','benchmark.tsv'),
  conda:
    '../envs/chimera.yml'
  shell:
    'python {params.script} --name chimera-1 --preprocess GraphST --adjacency SpaceFlow --train GraphST --clust mclust -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'

rule chimera2_clustering:
  input: 
    os.path.join('workflow','scripts','chimera_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','chimera-2','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','chimera_run.py'),
    n_domains = get_n_domains,
    outfile = os.path.join('results','{dataset}','chimera-2','{sample}','clustering.tsv'),
    outdir = os.path.join('results','{dataset}','chimera-2'),
  log:
    os.path.join('results','{dataset}','chimera-2','{sample}','chimera2.log'),
  benchmark:
    os.path.join('results','{dataset}','chimera-2','{sample}','benchmark.tsv'),
  conda:
    '../envs/chimera.yml'
  shell:
    'python {params.script} --name chimera-2 --preprocess GraphST --adjacency SCANIT --train GraphST --clust mclust -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'

rule chimera3_clustering:
  input: 
    os.path.join('workflow','scripts','chimera_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','chimera-3','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','chimera_run.py'),
    n_domains = get_n_domains,
    outfile = os.path.join('results','{dataset}','chimera-3','{sample}','clustering.tsv'),
    outdir = os.path.join('results','{dataset}','chimera-3'),
  log:
    os.path.join('results','{dataset}','chimera-3','{sample}','chimera3.log'),
  benchmark:
    os.path.join('results','{dataset}','chimera-3','{sample}','benchmark.tsv'),
  conda:
    '../envs/chimera.yml'
  shell:
    'python {params.script} --name chimera-3 --preprocess STAGATE --adjacency GraphST --train GraphST --clust mclust -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'

rule chimera4_clustering:
  input: 
    os.path.join('workflow','scripts','chimera_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','chimera-4','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','chimera_run.py'),
    n_domains = get_n_domains,
    outfile = os.path.join('results','{dataset}','chimera-4','{sample}','clustering.tsv'),
    outdir = os.path.join('results','{dataset}','chimera-4'),
  log:
    os.path.join('results','{dataset}','chimera-4','{sample}','chimera4.log'),
  benchmark:
    os.path.join('results','{dataset}','chimera-4','{sample}','benchmark.tsv'),
  conda:
    '../envs/chimera.yml'
  shell:
    'python {params.script} --name chimera-4 --preprocess STAGATE --adjacency CCST --train GraphST --clust mclust -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'

rule chimera5_clustering:
  input: 
    os.path.join('workflow','scripts','chimera_run.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','chimera-5','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','chimera_run.py'),
    n_domains = get_n_domains,
    outfile = os.path.join('results','{dataset}','chimera-5','{sample}','clustering.tsv'),
    outdir = os.path.join('results','{dataset}','chimera-5'),
  log:
    os.path.join('results','{dataset}','chimera-5','{sample}','chimera5.log'),
  benchmark:
    os.path.join('results','{dataset}','chimera-5','{sample}','benchmark.tsv'),
  conda:
    '../envs/chimera.yml'
  shell:
    'python {params.script} --name chimera-5 --preprocess GraphST --adjacency STAGATE --train GraphST --clust mclust -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
