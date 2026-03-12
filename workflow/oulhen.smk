import re
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

PCS = list(range(3, 13))
SAMPLES = ["integrated", "cont", "DAPT"]
PLOTFILES = ['ratio_group.png', 'Allstates.png', 'Freq_Prob_Energy.png', 'h.png', 'J.png', 'Basin.png', 'StatusNetwork_Subgraph.png', 'StatusNetwork_Subgraph_legend.png', 'StatusNetwork_Energy.png', 'StatusNetwork_Energy_legend.png', 'StatusNetwork_Ratio.png', 'StatusNetwork_Ratio_legend.png', 'StatusNetwork_State.png', 'StatusNetwork_State_legend.png', 'Landscape.png', 'discon_graph_1.png', 'discon_graph_2.png']
KFOLD = 5
SEED = 1
HPBASE_COMPARISONS = ['cont_DAPT', 'cont_DAPT_cov']

container: 'docker://koki/urchin_workflow_seurat:20251014'

rule all:
    input:
        # PC AUC Rank
        'output/oulhen/pc_auc_ovr.tsv',
        'output/oulhen/pc_auc_max.tsv',
        # PC-independent
        'plot/oulhen/integrated/pca_landscaper.png',
        expand('plot/oulhen/{sample}/dimplot_celltype_landscaper.png',
            sample=SAMPLES),
        'plot/oulhen/integrated/featureplot_ncount_rna_landscaper.png',
        expand('plot/oulhen/{sample}/cellular_density.png',
            sample=SAMPLES),
        # PC-dependent
        expand('plot/oulhen/{sample}/{pc}/sce.RData',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/energy.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/energy_hex.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/energy_rescaled.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/energy_rescaled_hex.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/energy_contour.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/cont_DAPT/{pc}/energy_diff.png',
            pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/basin.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/basin_piechart.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/states.png',
            sample=SAMPLES, pc=PCS),
        # Random Walk
        expand('plot/oulhen/{sample}/{pc}/P_metropolis_rw_Aboral_ectoderm.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_metropolis_rw_Oral_ectoderm.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_metropolis_rw_Ciliary_band.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_metropolis_rw_Neurons.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_glauber_rw_Aboral_ectoderm.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_glauber_rw_Oral_ectoderm.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_glauber_rw_Ciliary_band.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_glauber_rw_Neurons.png',
            sample=SAMPLES, pc=PCS),
        # Discon Graph
        expand('plot/oulhen/{sample}/{pc}/Landscaper/plot/discon_graph_1_new.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/Landscaper/plot/discon_graph_2_new.png',
            sample=SAMPLES, pc=PCS),
        # Graph Embedding
        expand('plot/oulhen/{sample}/{pc}/P_metropolis_emded.png',
            sample=SAMPLES, pc=PCS),
        expand('plot/oulhen/{sample}/{pc}/P_glauber_emded.png',
            sample=SAMPLES, pc=PCS),
        # Postprocess Landscaper
        expand('output/oulhen/{sample}/{pc}/seurat_landscaper.RData',
            sample=SAMPLES, pc=PCS),
        # Differential Landscape-associated Genes (DLGs)
        expand('output/oulhen/cont_DAPT/dlgs/{pc}/dlgs.tsv',
            pc=PCS),
        # DLGs Overlap (Oulhen vs HPBase)
        expand('output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}.tsv',
            pc=PCS, comparison=HPBASE_COMPARISONS),
        expand('plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}.png',
            pc=PCS, comparison=HPBASE_COMPARISONS),
        expand('plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}_venn_top.png',
            pc=PCS, comparison=HPBASE_COMPARISONS),
        expand('plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}_venn_bottom.png',
            pc=PCS, comparison=HPBASE_COMPARISONS),

#######################################
# Dimension Reduction & Preprocessing
#######################################
rule preprocess_landscaper_Oulhen:
    input:
        'data/Oulhen/Combined_Dim20_res0.5.rds'
    output:
        'output/oulhen/integrated/seurat_annotated_landscaper.RData',
        'output/oulhen/integrated/seurat.tsv',
        'output/oulhen/integrated/group.tsv',
        'output/oulhen/cont/group.tsv',
        'output/oulhen/DAPT/group.tsv'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_landscaper_Oulhen.txt'
    log:
        'logs/preprocess_landscaper_Oulhen.log'
    shell:
        'src/preprocess_landscaper_Oulhen.sh >& {log}'

rule stratify_seurat_Oulhen:
    input:
        'output/oulhen/integrated/seurat_annotated_landscaper.RData'
    output:
        'output/oulhen/cont/seurat_annotated_landscaper.RData',
        'output/oulhen/DAPT/seurat_annotated_landscaper.RData'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/stratify_seurat_Oulhen.txt'
    log:
        'logs/stratify_seurat_Oulhen.log'
    shell:
        'src/stratify_seurat_Oulhen.sh {input} {output} >& {log}'

rule preprocess_seurat_pcs:
    input:
        'output/oulhen/{sample}/seurat_annotated_landscaper.RData'
    output:
        'output/oulhen/{sample}/{pc}/seurat.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_seurat_pcs_oulhen_{sample}_{pc}.txt'
    log:
        'logs/preprocess_seurat_pcs_oulhen_{sample}_{pc}.log'
    shell:
        'src/preprocess_seurat_pcs.sh {input} {output} {wildcards.pc} >& {log}'

rule preprocess_binpca:
    input:
        'output/oulhen/{sample}/{pc}/seurat.tsv'
    output:
        'output/oulhen/{sample}/{pc}/binpca/BIN_DATA.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_binpca_oulhen_{sample}_{pc}.txt'
    log:
        'logs/preprocess_binpca_oulhen_{sample}_{pc}.log'
    shell:
        'src/preprocess_binpca.sh {input} {output} >& {log}'

rule pca_landscaper_Oulhen:
    input:
        'output/oulhen/integrated/seurat_annotated_landscaper.RData',
    output:
        'plot/oulhen/integrated/pca_landscaper.png'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/pca_landscaper_Oulhen.txt'
    log:
        'logs/pca_landscaper_Oulhen.log'
    shell:
        'src/pca_landscaper.sh {input} {output} >& {log}'

rule dimplot_celltype_landscaper_Oulhen:
    input:
        'output/oulhen/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/oulhen/{sample}/dimplot_celltype_landscaper.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_celltype_landscaper_Oulhen_{sample}.txt'
    log:
        'logs/dimplot_celltype_landscaper_Oulhen_{sample}.log'
    shell:
        'src/dimplot_celltype_landscaper_Oulhen.sh {input} {output} >& {log}'

rule featureplot_ncount_rna_landscaper_Oulhen:
    input:
        'output/oulhen/integrated/seurat_annotated_landscaper.RData'
    output:
        'plot/oulhen/integrated/featureplot_ncount_rna_landscaper.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_ncount_rna_landscaper_Oulhen.txt'
    log:
        'logs/featureplot_ncount_rna_landscaper_Oulhen.log'
    shell:
        'src/featureplot_ncount_rna_landscaper_Oulhen.sh {input} {output} >& {log}'

#######################################
# Landscaper
#######################################
rule landscaper_integrated:
    input:
        'output/oulhen/integrated/{pc}/binpca/BIN_DATA.tsv',
        'output/oulhen/integrated/group.tsv'
    output:
        'plot/oulhen/integrated/{pc}/Landscaper/plot/ratio_group.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/Allstates.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/h.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/J.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/Basin.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/Landscape.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/oulhen/integrated/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/oulhen/integrated/{pc}/Landscaper/Coordinate.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/Allstates.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/E.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/Basin.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/major_group.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/igraph.RData',
        'plot/oulhen/integrated/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/BIN_DATA'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_oulhen_integrated_{pc}.txt'
    log:
        'logs/landscaper_oulhen_integrated_{pc}.log'
    shell:
        'src/landscaper_integrated.sh {input} {output} >& {log}'

rule landscaper_cont:
    input:
        'output/oulhen/cont/{pc}/binpca/BIN_DATA.tsv',
        'output/oulhen/cont/group.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/Coordinate.tsv'
    output:
        'plot/oulhen/cont/{pc}/Landscaper/plot/ratio_group.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/Allstates.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/h.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/J.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/Basin.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/Landscape.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/oulhen/cont/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/oulhen/cont/{pc}/Landscaper/Basin.tsv',
        'plot/oulhen/cont/{pc}/Landscaper/major_group.tsv',
        'plot/oulhen/cont/{pc}/Landscaper/Allstates.tsv',
        'plot/oulhen/cont/{pc}/Landscaper/E.tsv',
        'plot/oulhen/cont/{pc}/Landscaper/igraph.RData',
        'plot/oulhen/cont/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/oulhen/cont/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/oulhen/cont/{pc}/Landscaper/BIN_DATA'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_oulhen_cont_{pc}.txt'
    log:
        'logs/landscaper_oulhen_cont_{pc}.log'
    shell:
        'src/landscaper_cont.sh {input} {output} >& {log}'

rule landscaper_DAPT:
    input:
        'output/oulhen/DAPT/{pc}/binpca/BIN_DATA.tsv',
        'output/oulhen/DAPT/group.tsv',
        'plot/oulhen/integrated/{pc}/Landscaper/Coordinate.tsv'
    output:
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/ratio_group.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/Allstates.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/h.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/J.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/Basin.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/Landscape.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/oulhen/DAPT/{pc}/Landscaper/Basin.tsv',
        'plot/oulhen/DAPT/{pc}/Landscaper/major_group.tsv',
        'plot/oulhen/DAPT/{pc}/Landscaper/Allstates.tsv',
        'plot/oulhen/DAPT/{pc}/Landscaper/E.tsv',
        'plot/oulhen/DAPT/{pc}/Landscaper/igraph.RData',
        'plot/oulhen/DAPT/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/oulhen/DAPT/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/oulhen/DAPT/{pc}/Landscaper/BIN_DATA'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_oulhen_DAPT_{pc}.txt'
    log:
        'logs/landscaper_oulhen_DAPT_{pc}.log'
    shell:
        'src/landscaper_DAPT.sh {input} {output} >& {log}'

#######################################
# Transition Probability
#######################################
rule transition_probability_Oulhen:
    input:
        'plot/oulhen/{sample}/{pc}/Landscaper/E.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/P_metropolis.tsv',
        'plot/oulhen/{sample}/{pc}/P_glauber.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/transition_probability_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/transition_probability_Oulhen_{sample}_{pc}.log'
    shell:
        'src/transition_probability_Oulhen.sh {input} {output} >& {log}'

#######################################
# Random Walk
#######################################
rule random_walk:
    input:
        'plot/oulhen/{sample}/{pc}/P_metropolis.tsv',
        'plot/oulhen/{sample}/{pc}/P_glauber.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/P_metropolis_rw.tsv',
        'plot/oulhen/{sample}/{pc}/P_glauber_rw.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/random_walk_oulhen_{sample}_{pc}.txt'
    log:
        'logs/random_walk_oulhen_{sample}_{pc}.log'
    shell:
        'src/random_walk.sh {input} {output} >& {log}'

rule plot_random_walk_metropolis_Oulhen:
    input:
        'plot/oulhen/integrated/{pc}/Landscaper/major_group.tsv',
        'plot/oulhen/{sample}/{pc}/P_metropolis_rw.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/P_metropolis_rw_Aboral_ectoderm.png',
        'plot/oulhen/{sample}/{pc}/P_metropolis_rw_Oral_ectoderm.png',
        'plot/oulhen/{sample}/{pc}/P_metropolis_rw_Ciliary_band.png',
        'plot/oulhen/{sample}/{pc}/P_metropolis_rw_Neurons.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_metropolis_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/plot_random_walk_metropolis_Oulhen_{sample}_{pc}.log'
    shell:
        'src/plot_random_walk_Oulhen.sh {input} {output} >& {log}'

rule plot_random_walk_glauber_Oulhen:
    input:
        'plot/oulhen/integrated/{pc}/Landscaper/major_group.tsv',
        'plot/oulhen/{sample}/{pc}/P_glauber_rw.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/P_glauber_rw_Aboral_ectoderm.png',
        'plot/oulhen/{sample}/{pc}/P_glauber_rw_Oral_ectoderm.png',
        'plot/oulhen/{sample}/{pc}/P_glauber_rw_Ciliary_band.png',
        'plot/oulhen/{sample}/{pc}/P_glauber_rw_Neurons.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_glauber_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/plot_random_walk_glauber_Oulhen_{sample}_{pc}.log'
    shell:
        'src/plot_random_walk_Oulhen.sh {input} {output} >& {log}'

#######################################
# Graph Embedding
#######################################
rule graph_embedding:
    input:
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates.tsv',
        'plot/oulhen/{sample}/{pc}/P_metropolis.tsv',
        'plot/oulhen/{sample}/{pc}/P_glauber.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/P_metropolis_emded.RData',
        'plot/oulhen/{sample}/{pc}/P_glauber_emded.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/graph_embedding_oulhen_{sample}_{pc}.txt'
    log:
        'logs/graph_embedding_oulhen_{sample}_{pc}.log'
    shell:
        'src/graph_embedding_Oulhen.sh {input} {output} >& {log}'

rule plot_graph_embedding:
    input:
        'plot/oulhen/{sample}/{pc}/Landscaper/E.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Basin.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/igraph.RData',
        'plot/oulhen/{sample}/{pc}/P_metropolis_emded.RData',
        'plot/oulhen/{sample}/{pc}/P_glauber_emded.RData'
    output:
        'plot/oulhen/{sample}/{pc}/P_metropolis_emded.png',
        'plot/oulhen/{sample}/{pc}/P_glauber_emded.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_graph_embedding_oulhen_{sample}_{pc}.txt'
    log:
        'logs/plot_graph_embedding_oulhen_{sample}_{pc}.log'
    shell:
        'src/plot_graph_embedding_Oulhen.sh {input} {output} >& {log}'

#######################################
# Tipping State & Dendrogram
#######################################
rule tipping_state:
    input:
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/E.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Basin.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/Landscaper/TippingState.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/tipping_state_oulhen_{sample}_{pc}.txt'
    log:
        'logs/tipping_state_oulhen_{sample}_{pc}.log'
    shell:
        'src/tipping_state_Oulhen.sh {input} {output} >& {log}'

rule dendrogram:
    input:
        'plot/oulhen/{sample}/{pc}/Landscaper/E.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Basin.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/TippingState.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/Landscaper/dendrogram_new.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dendrogram_oulhen_{sample}_{pc}.txt'
    log:
        'logs/dendrogram_oulhen_{sample}_{pc}.log'
    shell:
        'src/dendrogram_Oulhen.sh {input} {output} >& {log}'

rule plot_discon_graph:
    input:
        'plot/oulhen/{sample}/{pc}/Landscaper/dendrogram_new.RData'
    output:
        'plot/oulhen/{sample}/{pc}/Landscaper/plot/discon_graph_1_new.png',
        'plot/oulhen/{sample}/{pc}/Landscaper/plot/discon_graph_2_new.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_discon_graph_oulhen_{sample}_{pc}.txt'
    log:
        'logs/plot_discon_graph_oulhen_{sample}_{pc}.log'
    shell:
        'src/plot_discon_graph.sh {input} {output} >& {log}'

#######################################
# Featureplot Energy
#######################################
rule featureplot_energy_Oulhen:
    input:
        expand('plot/oulhen/{sample}/{{pc}}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/oulhen/{sample}/{pc}/sce.RData',
        'plot/oulhen/{sample}/{pc}/energy.png',
        'plot/oulhen/{sample}/{pc}/energy_hex.png',
        'plot/oulhen/{sample}/{pc}/energy_rescaled.png',
        'plot/oulhen/{sample}/{pc}/energy_rescaled_hex.png',
        'plot/oulhen/{sample}/{pc}/energy_contour.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/featureplot_energy_Oulhen_{sample}_{pc}.log'
    shell:
        'src/featureplot_energy_Oulhen_{wildcards.sample}.sh {wildcards.pc} {output} >& {log}'

rule featureplot_energy_Oulhen_cont_dapt:
    input:
        'plot/oulhen/integrated/{pc}/sce.RData',
        expand('plot/oulhen/cont/{{pc}}/Landscaper/plot/{p}', p=PLOTFILES),
        expand('plot/oulhen/DAPT/{{pc}}/Landscaper/plot/{p}', p=PLOTFILES)
    output:
        'plot/oulhen/cont_DAPT/{pc}/energy_diff.png'
    wildcard_constraints:
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_Oulhen_cont_DAPT_{pc}.txt'
    log:
        'logs/featureplot_energy_Oulhen_cont_DAPT_{pc}.log'
    shell:
        'src/featureplot_energy_Oulhen_cont_DAPT.sh {wildcards.pc} >& {log}'

#######################################
# Cellular Density
#######################################
rule plot_cellular_density:
    input:
        'output/oulhen/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/oulhen/{sample}/cellular_density.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_cellular_density_oulhen_{sample}.txt'
    log:
        'logs/plot_cellular_density_oulhen_{sample}.log'
    shell:
        'src/plot_cellular_density.sh {input} {output} >& {log}'

#######################################
# Dimplot Basin
#######################################
rule dimplot_basin_Oulhen:
    input:
        expand('plot/oulhen/{sample}/{{pc}}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/oulhen/{sample}/{pc}/basin.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_basin_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/dimplot_basin_Oulhen_{sample}_{pc}.log'
    shell:
        'src/dimplot_basin_Oulhen_{wildcards.sample}.sh {wildcards.pc} {output} >& {log}'

#######################################
# Basin Piechart
#######################################
rule piechart_basin_Oulhen:
    input:
        'output/oulhen/{sample}/seurat_annotated_landscaper.RData',
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/BIN_DATA',
        'plot/oulhen/{sample}/{pc}/Landscaper/Basin.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/basin_piechart.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/piechart_basin_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/piechart_basin_Oulhen_{sample}_{pc}.log'
    shell:
        'src/piechart_basin_Oulhen.sh {input} {output} >& {log}'

#######################################
# Dimplot States
#######################################
rule dimplot_states_Oulhen:
    input:
        expand('plot/oulhen/{sample}/{{pc}}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/oulhen/{sample}/{pc}/states.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_states_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/dimplot_states_Oulhen_{sample}_{pc}.log'
    shell:
        'src/dimplot_states_Oulhen_{wildcards.sample}.sh {wildcards.pc} {output} >& {log}'

#######################################
# Coarse-grained Vector Fields
#######################################
rule coarse_grained_vector_fields:
    input:
        'output/oulhen/integrated/seurat_annotated_landscaper.RData',
        'output/oulhen/{sample}/seurat_annotated_landscaper.RData',
        'output/oulhen/{sample}/{pc}/binpca/BIN_DATA.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates.tsv',
        'plot/oulhen/{sample}/{pc}/P_metropolis.tsv',
        'plot/oulhen/{sample}/{pc}/P_glauber.tsv'
    output:
        'plot/oulhen/{sample}/{pc}/P_metropolis_cg.RData',
        'plot/oulhen/{sample}/{pc}/P_glauber_cg.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/coarse_grained_vector_fields_oulhen_{sample}_{pc}.txt'
    log:
        'logs/coarse_grained_vector_fields_oulhen_{sample}_{pc}.log'
    shell:
        'src/coarse_grained_vector_fields_Oulhen.sh {input} {output} >& {log}'

#######################################
# PC AUC Rank
#######################################
rule pc_auc_rank:
    input:
        pca = 'output/oulhen/integrated/seurat.tsv',
        group = 'output/oulhen/integrated/group.tsv'
    output:
        ovr = 'output/oulhen/pc_auc_ovr.tsv',
        max = 'output/oulhen/pc_auc_max.tsv'
    params:
        kfold = KFOLD,
        seed = SEED
    resources:
        mem_mb = 10000000
    benchmark:
        'benchmarks/pc_auc_rank_oulhen.txt'
    log:
        'logs/pc_auc_rank_oulhen.log'
    shell:
        'src/pc_auc_rank.sh {input.pca} {input.group} {output.ovr} {output.max} {params.kfold} {params.seed} >& {log}'

#######################################
# Postprocess Landscaper
#######################################
rule postprocess_landscaper_Oulhen:
    input:
        'output/oulhen/{sample}/seurat_annotated_landscaper.RData',
        'output/oulhen/{sample}/{pc}/binpca/BIN_DATA.tsv',
        'output/oulhen/{sample}/group.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/E.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Basin.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/major_group.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Coordinate.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/igraph.RData',
        'plot/oulhen/{sample}/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/oulhen/{sample}/{pc}/Landscaper/Allstates_major_group.tsv'
    output:
        'output/oulhen/{sample}/{pc}/seurat_landscaper.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        pc='\\d+'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/postprocess_landscaper_Oulhen_{sample}_{pc}.txt'
    log:
        'logs/postprocess_landscaper_Oulhen_{sample}_{pc}.log'
    shell:
        'src/postprocess_landscaper_Oulhen.sh {input} {output} >& {log}'

#######################################
# Differential Landscape-associated Genes (DLGs)
#######################################
rule dlgs_cont_DAPT_Oulhen:
    input:
        'output/oulhen/cont/{pc}/seurat_landscaper.RData',
        'output/oulhen/DAPT/{pc}/seurat_landscaper.RData'
    output:
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs.tsv'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dlgs_cont_DAPT_Oulhen_{pc}.txt'
    log:
        'logs/dlgs_cont_DAPT_Oulhen_{pc}.log'
    shell:
        'src/dlgs.sh {input} {output} {wildcards.pc} >& {log}'

#######################################
# DLGs Overlap (Oulhen vs HPBase)
#######################################
rule dlgs_overlap:
    input:
        oulhen='output/oulhen/cont_DAPT/dlgs/{pc}/dlgs.tsv',
        hpbase='output/{comparison}/dlgs/{pc}/dlgs.tsv',
        annot='data/hpbase/HpulGenome_v1_annot.xlsx'
    output:
        tsv='output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}.tsv',
        png='plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}.png'
    wildcard_constraints:
        comparison='cont_DAPT|cont_DAPT_cov',
        pc='\\d+'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_overlap_{comparison}_{pc}.txt'
    log:
        'logs/dlgs_overlap_{comparison}_{pc}.log'
    shell:
        'src/dlgs_overlap.sh {input.oulhen} {input.hpbase} {input.annot} {output.tsv} {output.png} >& {log}'

rule dlgs_overlap_venn:
    input:
        overlap='output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}.tsv',
        oulhen='output/oulhen/cont_DAPT/dlgs/{pc}/dlgs.tsv',
        hpbase='output/{comparison}/dlgs/{pc}/dlgs.tsv'
    output:
        top='plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}_venn_top.png',
        bottom='plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_overlap_{comparison}_venn_bottom.png'
    wildcard_constraints:
        comparison='cont_DAPT|cont_DAPT_cov',
        pc='\\d+'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_overlap_venn_{comparison}_{pc}.txt'
    log:
        'logs/dlgs_overlap_venn_{comparison}_{pc}.log'
    shell:
        'src/dlgs_overlap_venn.sh {input.overlap} {input.oulhen} {input.hpbase} {output.top} {output.bottom} >& {log}'
