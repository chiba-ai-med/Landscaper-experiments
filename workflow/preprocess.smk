import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

PCS = list(range(2, 13))
SAMPLES = ['cont', 'DAPT', 'integrated', 'cont_cov', 'DAPT_cov', 'integrated_cov']

container: 'docker://koki/urchin_workflow_seurat:20251014'

rule all:
    input:
        # STRATIFY_DATA
        expand('output/{sample}/{pc}/binpca/BIN_DATA.tsv', sample=SAMPLES, pc=PCS),
        # PLOT
        expand('plot/{sample}/pca_landscaper.png', sample=SAMPLES),        
        expand('plot/{sample}/featureplot_ncount_rna_landscaper.png', sample=SAMPLES),
        expand('plot/{sample}/dimplot_celltype_landscaper.png', sample=SAMPLES),
        expand('plot/{sample}/dimplot_time_landscaper.png', sample=SAMPLES),
        expand('plot/{sample}/bindata/FINISH', sample=SAMPLES)

#######################################
# Dimension Reduction & Binarization
#######################################
rule preprocess_landscaper:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        # Stratified Data
        'output/integrated/seurat_annotated_landscaper.RData',
        'output/integrated_cov/seurat_annotated_landscaper.RData',
        'output/cont/seurat_annotated_landscaper.RData',
        'output/cont_cov/seurat_annotated_landscaper.RData',
        'output/DAPT/seurat_annotated_landscaper.RData',
        'output/DAPT_cov/seurat_annotated_landscaper.RData',
        # Group
        'output/integrated/group.tsv',
        'output/integrated_cov/group.tsv',
        'output/cont/group.tsv',
        'output/cont_cov/group.tsv',
        'output/DAPT/group.tsv',
        'output/DAPT_cov/group.tsv',
        # Covariates
        'output/integrated/cov.tsv',
        'output/integrated_cov/cov.tsv',
        'output/cont/cov.tsv',
        'output/cont_cov/cov.tsv',
        'output/DAPT/cov.tsv',
        'output/DAPT_cov/cov.tsv'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_landscaper.txt'
    log:
        'logs/preprocess_landscaper.log'
    shell:
        'src/preprocess_landscaper.sh >& {log}'

rule preprocess_seurat_pcs:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData'
    output:
        'output/{sample}/{pc}/seurat.tsv'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_seurat_pcs_{sample}_{pc}.txt'
    log:
        'logs/preprocess_seurat_pcs_{sample}_{pc}.log'
    shell:
        'src/preprocess_seurat_pcs.sh {input} {output} {wildcards.pc} >& {log}'

rule preprocess_binpca:
    input:
        'output/{sample}/{pc}/seurat.tsv'
    output:
        'output/{sample}/{pc}/binpca/BIN_DATA.tsv'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_binpca_{sample}_{pc}.txt'
    log:
        'logs/preprocess_binpca_{sample}_{pc}.log'
    shell:
        'src/preprocess_binpca.sh {input} {output} >& {log}'

#######################################
# Plot
#######################################
rule pca_landscaper:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/{sample}/pca_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/pca_landscaper_{sample}.txt'
    log:
        'logs/pca_landscaper_{sample}.log'
    shell:
        'src/pca_landscaper.sh {input} {output} >& {log}'

rule featureplot_ncount_rna_landscaper:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/{sample}/featureplot_ncount_rna_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_ncount_rna_landscaper_{sample}.txt'
    log:
        'logs/featureplot_ncount_rna_landscaper_{sample}.log'
    shell:
        'src/featureplot_ncount_rna_landscaper.sh {input} {output} >& {log}'

rule dimplot_celltype_landscaper:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/{sample}/dimplot_celltype_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_celltype_landscaper_{sample}.txt'
    log:
        'logs/dimplot_celltype_landscaper_{sample}.log'
    shell:
        'src/dimplot_celltype_landscaper.sh {input} {output} >& {log}'

rule dimplot_time_landscaper:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/{sample}/dimplot_time_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_time_landscaper_{sample}.txt'
    log:
        'logs/dimplot_time_landscaper_{sample}.log'
    shell:
        'src/dimplot_time_landscaper.sh {input} {output} >& {log}'

rule dimplot_bindata:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData',
        'output/{sample}/12/binpca/BIN_DATA.tsv'
    output:
        'plot/{sample}/bindata/FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_bindata_{sample}.txt'
    log:
        'logs/dimplot_bindata_{sample}.log'
    shell:
        'src/dimplot_bindata_{wildcards.sample}.sh {input} {output} >& {log}'
