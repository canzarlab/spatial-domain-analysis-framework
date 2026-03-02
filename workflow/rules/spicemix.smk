rule spicemix_preprocess:
  input:
    os.path.join('workflow','scripts','spicemix_prep.py'),
    os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
    os.path.join('resources','{dataset}','{sample}_coords.csv'),
  output:
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','expression_{sample}.txt'),
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','neighborhood_{sample}.txt'),
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','meta_{sample}.csv'),
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','genes_{sample}.txt'),
  params:
    preprocessing_script = os.path.join('workflow','scripts','spicemix_prep.py'),
    n_domains = get_n_domains,
  log:
    os.path.join('results','{dataset}','SpiceMix','{sample}','prep.log')
  conda:
    '../envs/spicemix_processing.yml'
  shell:
    'python {params.preprocessing_script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'


rule spicemix_clustering:
  input: 
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','expression_{sample}.txt'),
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','neighborhood_{sample}.txt'),
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','meta_{sample}.csv'),
    os.path.join('results','{dataset}','SpiceMix','{sample}','files','genes_{sample}.txt'),
    os.path.join('workflow','scripts','spicemix_run.py'),
  output:
    os.path.join('results','{dataset}','SpiceMix','{sample}','clustering.tsv'),
  params:
    script = os.path.join('workflow','scripts','spicemix_run.py'),
    n_domains = get_n_domains,
  log: 
    os.path.join('results','{dataset}','SpiceMix','{sample}','spicemix.log')
  benchmark:
    os.path.join('results','{dataset}','SpiceMix','{sample}','benchmark.tsv')
  conda:
    os.path.join('..','envs','spicemix.yml') 
  priority: -1
  shell:
    'python {params.script} -d {wildcards.dataset} -s {wildcards.sample} -n {params.n_domains}'
    '&> {log}'
#    '|| true ; touch {output}'
