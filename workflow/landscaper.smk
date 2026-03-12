import pandas as pd
import re
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

PCS = list(range(3, 13))
SAMPLES = ['cont', 'DAPT', 'integrated', 'cont_cov', 'DAPT_cov', 'integrated_cov']
SAMPLES2 = ['cont', 'DAPT', 'cont_cov', 'DAPT_cov']
GOS = ['bp', 'mf', 'cc']
LFILES = ['plot/ratio_group.png', 'plot/Allstates.png', 'plot/Freq_Prob_Energy.png',
    'plot/h.png', 'plot/J.png', 'plot/Basin.png',
    'plot/StatusNetwork_Subgraph.png', 'plot/StatusNetwork_Subgraph_legend.png',
    'plot/StatusNetwork_Energy.png', 'plot/StatusNetwork_Energy_legend.png',
    'plot/StatusNetwork_Ratio.png', 'plot/StatusNetwork_Ratio_legend.png',
    'plot/StatusNetwork_State.png', 'plot/StatusNetwork_State_legend.png',
    'plot/Landscape.png', 'plot/discon_graph_1.png', 'plot/discon_graph_2.png',
    'Coordinate.tsv', 'Allstates.tsv', 'E.tsv', 'Basin.tsv', 'major_group.tsv',
    'igraph.RData', 'EnergyBarrier.tsv', 'Allstates_major_group.tsv']
PLOTFILES = [f[5:] for f in LFILES if f.startswith('plot/')]

container: 'docker://koki/urchin_workflow_seurat:20251014'

rule all:
    input:
        # Landscaper
        expand('plot/{sample}/{pc}/Landscaper/{lfile}', sample=SAMPLES, pc=PCS, lfile=LFILES),
        # Postprocess
        expand('output/{sample}/{pc}/seurat_landscaper.RData', sample=SAMPLES, pc=PCS),
        # Downstream (per-sample, no PC)
        expand('plot/{sample}/cellular_density.png', sample=SAMPLES),
        # Downstream (per-sample, per-PC)
        expand('plot/{sample}/{pc}/energy.png', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/basin.png', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/basin_piechart.png', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/basin_trendtest_glm.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/states.png', sample=SAMPLES, pc=PCS),
        # Downstream (cross-sample, per-PC)
        expand('plot/integrated/{pc}/h.png', pc=PCS),
        expand('plot/integrated/{pc}/J.png', pc=PCS),
        expand('plot/cont_DAPT/{pc}/energy_diff.png', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/energy_diff.png', pc=PCS),
        expand('plot/{sample2}/{pc}/landscape_rescaled.png', sample2=SAMPLES2, pc=PCS),
        expand('plot/{pc}/stable_states.png', pc=PCS),
        # Differential Landscape-associated Genes (DLGs)
        expand('output/cont_DAPT/dlgs/{pc}/dlgs.tsv',
            pc=PCS),
        expand('output/cont_DAPT_cov/dlgs/{pc}/dlgs.tsv',
            pc=PCS)

#######################################
# Energy Landscape Analysis
#######################################
rule landscaper_integrated:
    input:
        'output/integrated/{pc}/binpca/BIN_DATA.tsv',
        'output/integrated/group.tsv'
    output:
        'plot/integrated/{pc}/Landscaper/plot/ratio_group.png',
        'plot/integrated/{pc}/Landscaper/plot/Allstates.png',
        'plot/integrated/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/integrated/{pc}/Landscaper/plot/h.png',
        'plot/integrated/{pc}/Landscaper/plot/J.png',
        'plot/integrated/{pc}/Landscaper/plot/Basin.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/integrated/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/integrated/{pc}/Landscaper/plot/Landscape.png',
        'plot/integrated/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/integrated/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/integrated/{pc}/Landscaper/Coordinate.tsv',
        'plot/integrated/{pc}/Landscaper/Allstates.tsv',
        'plot/integrated/{pc}/Landscaper/E.tsv',
        'plot/integrated/{pc}/Landscaper/Basin.tsv',
        'plot/integrated/{pc}/Landscaper/major_group.tsv',
        'plot/integrated/{pc}/Landscaper/igraph.RData',
        'plot/integrated/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/integrated/{pc}/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_integrated_{pc}.txt'
    log:
        'logs/landscaper_integrated_{pc}.log'
    shell:
        'src/landscaper_integrated.sh {input} {output} >& {log}'

rule landscaper_cont:
    input:
        'output/cont/{pc}/binpca/BIN_DATA.tsv',
        'output/cont/group.tsv',
        'plot/integrated/{pc}/Landscaper/Coordinate.tsv'
    output:
        'plot/cont/{pc}/Landscaper/plot/ratio_group.png',
        'plot/cont/{pc}/Landscaper/plot/Allstates.png',
        'plot/cont/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/cont/{pc}/Landscaper/plot/h.png',
        'plot/cont/{pc}/Landscaper/plot/J.png',
        'plot/cont/{pc}/Landscaper/plot/Basin.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/cont/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/cont/{pc}/Landscaper/plot/Landscape.png',
        'plot/cont/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/cont/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/cont/{pc}/Landscaper/Basin.tsv',
        'plot/cont/{pc}/Landscaper/major_group.tsv',
        'plot/cont/{pc}/Landscaper/Coordinate.tsv',
        'plot/cont/{pc}/Landscaper/Allstates.tsv',
        'plot/cont/{pc}/Landscaper/E.tsv',
        'plot/cont/{pc}/Landscaper/igraph.RData',
        'plot/cont/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/cont/{pc}/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_cont_{pc}.txt'
    log:
        'logs/landscaper_cont_{pc}.log'
    shell:
        'src/landscaper_cont.sh {input} {output} >& {log}'

rule landscaper_DAPT:
    input:
        'output/DAPT/{pc}/binpca/BIN_DATA.tsv',
        'output/DAPT/group.tsv',
        'plot/integrated/{pc}/Landscaper/Coordinate.tsv'
    output:
        'plot/DAPT/{pc}/Landscaper/plot/ratio_group.png',
        'plot/DAPT/{pc}/Landscaper/plot/Allstates.png',
        'plot/DAPT/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/DAPT/{pc}/Landscaper/plot/h.png',
        'plot/DAPT/{pc}/Landscaper/plot/J.png',
        'plot/DAPT/{pc}/Landscaper/plot/Basin.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/DAPT/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/DAPT/{pc}/Landscaper/plot/Landscape.png',
        'plot/DAPT/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/DAPT/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/DAPT/{pc}/Landscaper/Basin.tsv',
        'plot/DAPT/{pc}/Landscaper/major_group.tsv',
        'plot/DAPT/{pc}/Landscaper/Coordinate.tsv',
        'plot/DAPT/{pc}/Landscaper/Allstates.tsv',
        'plot/DAPT/{pc}/Landscaper/E.tsv',
        'plot/DAPT/{pc}/Landscaper/igraph.RData',
        'plot/DAPT/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/DAPT/{pc}/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_DAPT_{pc}.txt'
    log:
        'logs/landscaper_DAPT_{pc}.log'
    shell:
        'src/landscaper_DAPT.sh {input} {output} >& {log}'

rule landscaper_integrated_cov:
    input:
        'output/integrated_cov/{pc}/binpca/BIN_DATA.tsv',
        'output/integrated_cov/group.tsv',
        'plot/integrated/{pc}/Landscaper/Coordinate.tsv',
        'output/integrated_cov/cov.tsv'
    output:
        'plot/integrated_cov/{pc}/Landscaper/plot/ratio_group.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/Allstates.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/h.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/J.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/Basin.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/Landscape.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/integrated_cov/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/integrated_cov/{pc}/Landscaper/Coordinate.tsv',
        'plot/integrated_cov/{pc}/Landscaper/Allstates.tsv',
        'plot/integrated_cov/{pc}/Landscaper/E.tsv',
        'plot/integrated_cov/{pc}/Landscaper/Basin.tsv',
        'plot/integrated_cov/{pc}/Landscaper/major_group.tsv',
        'plot/integrated_cov/{pc}/Landscaper/igraph.RData',
        'plot/integrated_cov/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/integrated_cov/{pc}/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_integrated_cov_{pc}.txt'
    log:
        'logs/landscaper_integrated_cov_{pc}.log'
    shell:
        'src/landscaper_integrated_cov.sh {input} {output} >& {log}'

rule landscaper_cont_cov:
    input:
        'output/cont_cov/{pc}/binpca/BIN_DATA.tsv',
        'output/cont_cov/group.tsv',
        'plot/integrated/{pc}/Landscaper/Coordinate.tsv',
        'output/cont_cov/cov.tsv'
    output:
        'plot/cont_cov/{pc}/Landscaper/plot/ratio_group.png',
        'plot/cont_cov/{pc}/Landscaper/plot/Allstates.png',
        'plot/cont_cov/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/cont_cov/{pc}/Landscaper/plot/h.png',
        'plot/cont_cov/{pc}/Landscaper/plot/J.png',
        'plot/cont_cov/{pc}/Landscaper/plot/Basin.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/cont_cov/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/cont_cov/{pc}/Landscaper/plot/Landscape.png',
        'plot/cont_cov/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/cont_cov/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/cont_cov/{pc}/Landscaper/h.tsv',
        'plot/cont_cov/{pc}/Landscaper/J.tsv',
        'plot/cont_cov/{pc}/Landscaper/g.txt',
        'plot/cont_cov/{pc}/Landscaper/Basin.tsv',
        'plot/cont_cov/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/cont_cov/{pc}/Landscaper/major_group.tsv',
        'plot/cont_cov/{pc}/Landscaper/Coordinate.tsv',
        'plot/cont_cov/{pc}/Landscaper/Allstates.tsv',
        'plot/cont_cov/{pc}/Landscaper/E.tsv',
        'plot/cont_cov/{pc}/Landscaper/igraph.RData',
        'plot/cont_cov/{pc}/Landscaper/EnergyBarrier.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_cont_cov_{pc}.txt'
    log:
        'logs/landscaper_cont_cov_{pc}.log'
    shell:
        'src/landscaper_cont_cov.sh {input} {output} >& {log}'

rule landscaper_DAPT_cov:
    input:
        'output/DAPT_cov/{pc}/binpca/BIN_DATA.tsv',
        'output/DAPT_cov/group.tsv',
        'plot/integrated/{pc}/Landscaper/Coordinate.tsv',
        'output/DAPT_cov/cov.tsv'
    output:
        'plot/DAPT_cov/{pc}/Landscaper/plot/ratio_group.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/Allstates.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/h.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/J.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/Basin.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_State.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/Landscape.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/discon_graph_1.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/discon_graph_2.png',
        'plot/DAPT_cov/{pc}/Landscaper/h.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/J.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/g.txt',
        'plot/DAPT_cov/{pc}/Landscaper/Basin.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/major_group.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/Coordinate.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/Allstates.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/E.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/igraph.RData',
        'plot/DAPT_cov/{pc}/Landscaper/EnergyBarrier.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_DAPT_cov_{pc}.txt'
    log:
        'logs/landscaper_DAPT_cov_{pc}.log'
    shell:
        'src/landscaper_DAPT_cov.sh {input} {output} >& {log}'

#######################################
# Postprocess: consolidate into Seurat
#######################################
rule postprocess_landscaper:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData',
        'output/{sample}/{pc}/binpca/BIN_DATA.tsv',
        'output/{sample}/group.tsv',
        'output/{sample}/cov.tsv',
        'plot/{sample}/{pc}/Landscaper/Allstates.tsv',
        'plot/{sample}/{pc}/Landscaper/E.tsv',
        'plot/{sample}/{pc}/Landscaper/Basin.tsv',
        'plot/{sample}/{pc}/Landscaper/major_group.tsv',
        'plot/{sample}/{pc}/Landscaper/Coordinate.tsv',
        'plot/{sample}/{pc}/Landscaper/igraph.RData',
        'plot/{sample}/{pc}/Landscaper/EnergyBarrier.tsv',
        'plot/{sample}/{pc}/Landscaper/Allstates_major_group.tsv'
    output:
        'output/{sample}/{pc}/seurat_landscaper.RData'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/postprocess_landscaper_{sample}_{pc}.txt'
    log:
        'logs/postprocess_landscaper_{sample}_{pc}.log'
    shell:
        'src/postprocess_landscaper.sh {input} {output} >& {log}'

#######################################
# Downstream Analysis
#######################################
rule plot_h:
    input:
        'plot/cont/{pc}/Landscaper/plot/h.png',
        'plot/DAPT/{pc}/Landscaper/plot/h.png',
        'plot/cont_cov/{pc}/Landscaper/plot/h.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/h.png'
    output:
        'plot/integrated/{pc}/h.png',
        'plot/integrated_cov/{pc}/h.png'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_h_{pc}.txt'
    log:
        'logs/plot_h_{pc}.log'
    shell:
        'src/plot_h.sh {wildcards.pc} {output} >& {log}'

rule plot_J:
    input:
        'plot/cont/{pc}/Landscaper/plot/J.png',
        'plot/DAPT/{pc}/Landscaper/plot/J.png',
        'plot/cont_cov/{pc}/Landscaper/plot/J.png',
        'plot/DAPT_cov/{pc}/Landscaper/plot/J.png'
    output:
        'plot/integrated/{pc}/J.png',
        'plot/integrated_cov/{pc}/J.png'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_J_{pc}.txt'
    log:
        'logs/plot_J_{pc}.log'
    shell:
        'src/plot_J.sh {wildcards.pc} {output} >& {log}'

rule featureplot_energy:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/sce.RData',
        'plot/{sample}/{pc}/energy.png',
        'plot/{sample}/{pc}/energy_hex.png',
        'plot/{sample}/{pc}/energy_rescaled.png',
        'plot/{sample}/{pc}/energy_rescaled_hex.png',
        'plot/{sample}/{pc}/energy_splitby.png',
        'plot/{sample}/{pc}/energy_contour.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_{sample}_{pc}.txt'
    log:
        'logs/featureplot_energy_{sample}_{pc}.log'
    shell:
        'src/featureplot_energy.sh {input} {output} >& {log}'

rule featureplot_energy_cont_dapt:
    input:
        'plot/integrated/{pc}/sce.RData',
        'output/cont/{pc}/seurat_landscaper.RData',
        'output/DAPT/{pc}/seurat_landscaper.RData'
    output:
        'plot/cont_DAPT/{pc}/energy_diff.png'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_cont_DAPT_{pc}.txt'
    log:
        'logs/featureplot_energy_cont_DAPT_{pc}.log'
    shell:
        'src/featureplot_energy_cont_dapt.sh {input} {output} >& {log}'

rule featureplot_energy_cont_dapt_cov:
    input:
        'plot/integrated_cov/{pc}/sce.RData',
        'output/cont_cov/{pc}/seurat_landscaper.RData',
        'output/DAPT_cov/{pc}/seurat_landscaper.RData'
    output:
        'plot/cont_DAPT_cov/{pc}/energy_diff.png'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_cont_DAPT_cov_{pc}.txt'
    log:
        'logs/featureplot_energy_cont_DAPT_cov_{pc}.log'
    shell:
        'src/featureplot_energy_cont_dapt.sh {input} {output} >& {log}'

rule plot_cellular_density:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/{sample}/cellular_density.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_cellular_density_{sample}.txt'
    log:
        'logs/plot_cellular_density_{sample}.log'
    shell:
        'src/plot_cellular_density.sh {input} {output} >& {log}'

rule dimplot_basin:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/basin.png',
        'plot/{sample}/{pc}/basin_splitby.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_basin_{sample}_{pc}.txt'
    log:
        'logs/dimplot_basin_{sample}_{pc}.log'
    shell:
        'src/dimplot_basin.sh {input} {output} >& {log}'

rule piechart_basin:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData',
        'plot/{sample}/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/{sample}/{pc}/Landscaper/BIN_DATA',
        'plot/{sample}/{pc}/Landscaper/Basin.tsv'
    output:
        'plot/{sample}/{pc}/basin_piechart.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/piechart_basin_{sample}_{pc}.txt'
    log:
        'logs/piechart_basin_{sample}_{pc}.log'
    shell:
        'src/piechart_basin.sh {input} {output} >& {log}'

rule trendtest_basin:
    input:
        'output/{sample}/seurat_annotated_landscaper.RData',
        'plot/{sample}/{pc}/Landscaper/Allstates_major_group.tsv',
        'plot/{sample}/{pc}/Landscaper/BIN_DATA',
        'plot/{sample}/{pc}/Landscaper/Basin.tsv'
    output:
        'plot/{sample}/{pc}/basin_trendtest_glm.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/trendtest_basin_{sample}_{pc}.txt'
    log:
        'logs/trendtest_basin_{sample}_{pc}.log'
    shell:
        'src/trendtest_basin.sh {input} {output} >& {log}'

rule dimplot_states:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/states.png',
        'plot/{sample}/{pc}/states_splitby.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_states_{sample}_{pc}.txt'
    log:
        'logs/dimplot_states_{sample}_{pc}.log'
    shell:
        'src/dimplot_states.sh {input} {output} >& {log}'

rule plot_landscape:
    input:
        expand('plot/{sample}/{{pc}}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/cont/{pc}/landscape_rescaled.png',
        'plot/DAPT/{pc}/landscape_rescaled.png',
        'plot/cont_cov/{pc}/landscape_rescaled.png',
        'plot/DAPT_cov/{pc}/landscape_rescaled.png'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_landscape_{pc}.txt'
    log:
        'logs/plot_landscape_{pc}.log'
    shell:
        'src/plot_landscape.sh {wildcards.pc} {output} >& {log}'

rule plot_stable_states:
    input:
        'plot/cont/{pc}/Landscaper/Basin.tsv',
        'plot/DAPT/{pc}/Landscaper/Basin.tsv',
        'plot/cont_cov/{pc}/Landscaper/Basin.tsv',
        'plot/DAPT_cov/{pc}/Landscaper/Basin.tsv'
    output:
        'plot/{pc}/stable_states.png'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_stable_states_{pc}.txt'
    log:
        'logs/plot_stable_states_{pc}.log'
    shell:
        'src/plot_stable_states.sh {wildcards.pc} {output} >& {log}'


# Differential Landscape-associated Genes (DLGs)
rule dlgs_cont_DAPT:
    input:
        'output/cont/{pc}/seurat_landscaper.RData',
        'output/DAPT/{pc}/seurat_landscaper.RData'
    output:
        'output/cont_DAPT/dlgs/{pc}/dlgs.tsv'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dlgs_cont_DAPT_{pc}.txt'
    log:
        'logs/dlgs_cont_DAPT_{pc}.log'
    shell:
        'src/dlgs.sh {input} {output} {wildcards.pc} >& {log}'

rule dlgs_cont_DAPT_cov:
    input:
        'output/cont_cov/{pc}/seurat_landscaper.RData',
        'output/DAPT_cov/{pc}/seurat_landscaper.RData'
    output:
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs.tsv'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dlgs_cont_DAPT_cov_{pc}.txt'
    log:
        'logs/dlgs_cont_DAPT_cov_{pc}.log'
    shell:
        'src/dlgs.sh {input} {output} {wildcards.pc} >& {log}'
