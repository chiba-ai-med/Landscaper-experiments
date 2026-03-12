import re
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

PCS = list(range(3, 13))
SAMPLES = ['cont', 'DAPT', 'integrated', 'cont_cov', 'DAPT_cov', 'integrated_cov']

container: 'docker://koki/urchin_workflow_seurat:20251014'

rule all:
    input:
        expand('plot/{sample}/{pc}/P_metropolis_emded.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_emded.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_metropolis_emded.png', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_emded.png', sample=SAMPLES, pc=PCS)

#######################################
# Graph Embedding
#######################################
rule graph_embedding:
    input:
        'plot/{sample}/{pc}/P_metropolis.tsv',
        'plot/{sample}/{pc}/P_glauber.tsv'
    output:
        'plot/{sample}/{pc}/P_metropolis_emded.RData',
        'plot/{sample}/{pc}/P_glauber_emded.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/graph_embedding_{sample}_{pc}.txt'
    log:
        'logs/graph_embedding_{sample}_{pc}.log'
    shell:
        'src/graph_embedding.sh {input} {output} >& {log}'

#######################################
# Plot Graph Embedding
#######################################
rule plot_graph_embedding:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis_emded.RData',
        'plot/{sample}/{pc}/P_glauber_emded.RData'
    output:
        'plot/{sample}/{pc}/P_metropolis_emded.png',
        'plot/{sample}/{pc}/P_glauber_emded.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_graph_embedding_{sample}_{pc}.txt'
    log:
        'logs/plot_graph_embedding_{sample}_{pc}.log'
    shell:
        'src/plot_graph_embedding.sh {input} {output} >& {log}'
