rule seurat_clustering:
    input:
        os.path.join("workflow", "scripts", "seurat_run.R"),
        os.path.join('resources','{dataset}','{sample}_UMImatrix.csv'),
        os.path.join('resources','{dataset}','{sample}_coords.csv'),
    output:
        os.path.join("results","{dataset}","Seurat","{sample}","clustering.tsv"),
    message:
        ""
    params:
        seurat_script = os.path.join("workflow", "scripts", "seurat_run.R"),
        n_cluster = lambda wildcards: sp_df[(sp_df.dataset_id == wildcards.dataset) & (sp_df.sample_id == wildcards.sample)].n_cluster.values[0],
    log:
        os.path.join("results", "{dataset}", "Seurat", "{sample}", "seurat.log")
    benchmark:
        os.path.join('results','{dataset}','Seurat','{sample}','benchmark.tsv')
    conda:
        "../envs/seurat.yml"
    shell:
        'Rscript {params.seurat_script} {wildcards.dataset} {wildcards.sample} {params.n_cluster}'
        '&> {log}'
        '|| true ; touch {output}'
