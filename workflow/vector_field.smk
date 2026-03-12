import re
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

PCS = list(range(3, 13))
SAMPLES = ['cont', 'DAPT', 'integrated', 'cont_cov', 'DAPT_cov', 'integrated_cov']
CG_MODES = ['', '_germlayer', '_cluster', '_sample', '_celltype', '_state', '_energy']

container: 'docker://koki/urchin_workflow_seurat:20251014'

rule all:
    input:
        expand('plot/{sample}/{pc}/P_metropolis_cg.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_cg.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_metropolis_cg{mode}.png',
            sample=SAMPLES, pc=PCS, mode=CG_MODES),
        expand('plot/{sample}/{pc}/P_glauber_cg{mode}.png',
            sample=SAMPLES, pc=PCS, mode=CG_MODES),
        expand('plot/{sample}/{pc}/P_metropolis_cg{mode}_streamline.png',
            sample=SAMPLES, pc=PCS, mode=CG_MODES),
        expand('plot/{sample}/{pc}/P_glauber_cg{mode}_streamline.png',
            sample=SAMPLES, pc=PCS, mode=CG_MODES)

#######################################
# Coarse-Grained Vector Fields
#######################################
rule coarse_grained_vector_fields:
    input:
        'output/integrated/{pc}/seurat_landscaper.RData',
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis.tsv',
        'plot/{sample}/{pc}/P_glauber.tsv'
    output:
        'plot/{sample}/{pc}/P_metropolis_cg.RData',
        'plot/{sample}/{pc}/P_glauber_cg.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/coarse_grained_vector_fields_{sample}_{pc}.txt'
    log:
        'logs/coarse_grained_vector_fields_{sample}_{pc}.log'
    shell:
        'src/coarse_grained_vector_fields.sh {input} {output} >& {log}'

#######################################
# Plot Coarse-Grained Vector Fields
#######################################
rule plot_coarse_grained_vector_fields:
    input:
        'output/integrated/{pc}/seurat_landscaper.RData',
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis_cg.RData',
        'plot/{sample}/{pc}/P_glauber_cg.RData'
    output:
        'plot/{sample}/{pc}/P_metropolis_cg.png',
        'plot/{sample}/{pc}/P_metropolis_cg_germlayer.png',
        'plot/{sample}/{pc}/P_metropolis_cg_cluster.png',
        'plot/{sample}/{pc}/P_metropolis_cg_sample.png',
        'plot/{sample}/{pc}/P_metropolis_cg_celltype.png',
        'plot/{sample}/{pc}/P_metropolis_cg_state.png',
        'plot/{sample}/{pc}/P_metropolis_cg_energy.png',
        'plot/{sample}/{pc}/P_glauber_cg.png',
        'plot/{sample}/{pc}/P_glauber_cg_germlayer.png',
        'plot/{sample}/{pc}/P_glauber_cg_cluster.png',
        'plot/{sample}/{pc}/P_glauber_cg_sample.png',
        'plot/{sample}/{pc}/P_glauber_cg_celltype.png',
        'plot/{sample}/{pc}/P_glauber_cg_state.png',
        'plot/{sample}/{pc}/P_glauber_cg_energy.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_coarse_grained_vector_fields_{sample}_{pc}.txt'
    log:
        'logs/plot_coarse_grained_vector_fields_{sample}_{pc}.log'
    shell:
        'src/plot_coarse_grained_vector_fields.sh {input} {output} >& {log}'

#######################################
# Plot Coarse-Grained Vector Fields (Streamline)
#######################################
rule plot_coarse_grained_vector_fields_streamline:
    input:
        'output/integrated/{pc}/seurat_landscaper.RData',
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis_cg.RData',
        'plot/{sample}/{pc}/P_glauber_cg.RData'
    output:
        'plot/{sample}/{pc}/P_metropolis_cg_streamline.png',
        'plot/{sample}/{pc}/P_metropolis_cg_germlayer_streamline.png',
        'plot/{sample}/{pc}/P_metropolis_cg_cluster_streamline.png',
        'plot/{sample}/{pc}/P_metropolis_cg_sample_streamline.png',
        'plot/{sample}/{pc}/P_metropolis_cg_celltype_streamline.png',
        'plot/{sample}/{pc}/P_metropolis_cg_state_streamline.png',
        'plot/{sample}/{pc}/P_metropolis_cg_energy_streamline.png',
        'plot/{sample}/{pc}/P_glauber_cg_streamline.png',
        'plot/{sample}/{pc}/P_glauber_cg_germlayer_streamline.png',
        'plot/{sample}/{pc}/P_glauber_cg_cluster_streamline.png',
        'plot/{sample}/{pc}/P_glauber_cg_sample_streamline.png',
        'plot/{sample}/{pc}/P_glauber_cg_celltype_streamline.png',
        'plot/{sample}/{pc}/P_glauber_cg_state_streamline.png',
        'plot/{sample}/{pc}/P_glauber_cg_energy_streamline.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_coarse_grained_vector_fields_streamline_{sample}_{pc}.txt'
    log:
        'logs/plot_coarse_grained_vector_fields_streamline_{sample}_{pc}.log'
    shell:
        'src/plot_coarse_grained_vector_fields_streamline.sh {input} {output} >& {log}'
