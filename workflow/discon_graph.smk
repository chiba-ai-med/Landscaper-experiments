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
        expand('plot/{sample}/{pc}/TippingState.tsv', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/dendrogram_new.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/discon_graph_1_new.png', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/discon_graph_2_new.png', sample=SAMPLES, pc=PCS)

#######################################
# Tipping State
#######################################
rule tipping_state:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/TippingState.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/tipping_state_{sample}_{pc}.txt'
    log:
        'logs/tipping_state_{sample}_{pc}.log'
    shell:
        'src/tipping_state.sh {input} {output} >& {log}'

#######################################
# Dendrogram
#######################################
rule dendrogram:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/TippingState.tsv'
    output:
        'plot/{sample}/{pc}/dendrogram_new.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dendrogram_{sample}_{pc}.txt'
    log:
        'logs/dendrogram_{sample}_{pc}.log'
    shell:
        'src/dendrogram.sh {input} {output} >& {log}'

#######################################
# Plot Disconnectivity Graph
#######################################
rule plot_discon_graph:
    input:
        'plot/{sample}/{pc}/dendrogram_new.RData'
    output:
        'plot/{sample}/{pc}/discon_graph_1_new.png',
        'plot/{sample}/{pc}/discon_graph_2_new.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_discon_graph_{sample}_{pc}.txt'
    log:
        'logs/plot_discon_graph_{sample}_{pc}.log'
    shell:
        'src/plot_discon_graph.sh {input} {output} >& {log}'
