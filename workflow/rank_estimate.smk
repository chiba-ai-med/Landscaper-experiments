from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

KFOLD = 5
SEED = 1
PCA_INPUT = "output/integrated/12/seurat.tsv"
GROUP_INPUT = "output/integrated/group.tsv"

container: 'docker://koki/urchin_workflow_seurat:20251014'

rule all:
    input:
        "output/pc_auc_ovr.tsv",
        "output/pc_auc_max.tsv"

rule pc_auc_rank:
    input:
        pca = PCA_INPUT,
        group = GROUP_INPUT
    output:
        ovr = "output/pc_auc_ovr.tsv",
        max = "output/pc_auc_max.tsv"
    params:
        kfold = KFOLD,
        seed = SEED
    resources:
        mem_mb = 10000000
    log:
        "logs/pc_auc_rank.log"
    benchmark:
        "benchmarks/pc_auc_rank.txt"
    shell:
        "src/pc_auc_rank.sh {input.pca} {input.group} {output.ovr} {output.max} {params.kfold} {params.seed} >& {log}"
