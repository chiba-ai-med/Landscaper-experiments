import re
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

PCS = list(range(3, 13))
SAMPLES = ['cont', 'DAPT', 'integrated', 'cont_cov', 'DAPT_cov', 'integrated_cov']
OULHEN_SAMPLES = ['cont', 'DAPT', 'integrated']
GOS = ['bp', 'mf', 'cc']

container: 'docker://koki/urchin_workflow_gsea:20260113'

rule all:
    input:
        expand('plot/{sample}/{pc}/basin_gsea_{go}.RData',
            sample=SAMPLES, pc=PCS, go=GOS),
        expand('plot/{sample}/{pc}/basin_gsea_{go}_plots/FINISH',
            sample=SAMPLES, pc=PCS, go=GOS),
        expand('plot/{sample}/{pc}/basin_gsea_{go}_hexbin/FINISH',
            sample=SAMPLES, pc=PCS, go=GOS),
        expand('plot/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}_plots/FINISH',
            pc=PCS, go=GOS),
        expand('plot/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_{go}_plots/FINISH',
            pc=PCS, go=GOS),
        # DLGs waterfall plots
        expand('plot/cont_DAPT/dlgs/{pc}/dlgs_waterfall.png', pc=PCS),
        expand('plot/cont_DAPT_cov/dlgs/{pc}/dlgs_waterfall.png', pc=PCS),
        # DLGs GSEA dot plots
        expand('plot/cont_DAPT/dlgs/{pc}/dlgs_gsea_dotplot.png', pc=PCS),
        expand('plot/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_dotplot.png', pc=PCS),
        # Oulhen
        expand('plot/oulhen/{sample}/{pc}/basin_gsea_{go}.RData',
            sample=OULHEN_SAMPLES, pc=PCS, go=GOS),
        expand('plot/oulhen/{sample}/{pc}/basin_gsea_{go}_plots/FINISH',
            sample=OULHEN_SAMPLES, pc=PCS, go=GOS),
        expand('plot/oulhen/{sample}/{pc}/basin_gsea_{go}_hexbin/FINISH',
            sample=OULHEN_SAMPLES, pc=PCS, go=GOS),
        expand('output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}.RData',
            pc=PCS, go=GOS),
        expand('plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}_plots/FINISH',
            pc=PCS, go=GOS),
        # DLGs waterfall plots (Oulhen)
        expand('plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_waterfall.png', pc=PCS),
        # DLGs GSEA dot plots (Oulhen)
        expand('plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_dotplot.png', pc=PCS)

#######################################
# GSEA Basin
#######################################
rule gsea_basin:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'data/go_{go}_hpbase.RData'
    output:
        'plot/{sample}/{pc}/basin_gsea_{go}.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_basin_{sample}_{pc}_{go}.txt'
    log:
        'logs/gsea_basin_{sample}_{pc}_{go}.log'
    shell:
        'src/gsea_basin.sh {input} {output} >& {log}'

#######################################
# GSEA Basin Plots
#######################################
rule gsea_basin_plots:
    input:
        'plot/{sample}/{pc}/basin_gsea_{go}.RData'
    output:
        'plot/{sample}/{pc}/basin_gsea_{go}_plots/FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_basin_plots_{sample}_{pc}_{go}.txt'
    log:
        'logs/gsea_basin_plots_{sample}_{pc}_{go}.log'
    shell:
        'src/gsea_basin_plots.sh {input} {output} >& {log}'

#######################################
# GSEA for DLGs
#######################################
rule gsea_dlgs_cont_DAPT:
    input:
        'output/cont_DAPT/dlgs/{pc}/dlgs.tsv',
        'data/go_{go}_hpbase.RData'
    output:
        'output/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}.RData'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_dlgs_cont_DAPT_{pc}_{go}.txt'
    log:
        'logs/gsea_dlgs_cont_DAPT_{pc}_{go}.log'
    shell:
        'src/gsea_dlgs.sh {input} {output} >& {log}'

rule gsea_dlgs_cont_DAPT_cov:
    input:
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs.tsv',
        'data/go_{go}_hpbase.RData'
    output:
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_{go}.RData'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_dlgs_cont_DAPT_cov_{pc}_{go}.txt'
    log:
        'logs/gsea_dlgs_cont_DAPT_cov_{pc}_{go}.log'
    shell:
        'src/gsea_dlgs.sh {input} {output} >& {log}'

#######################################
# GSEA Plots for DLGs
#######################################
rule gsea_dlgs_plots_cont_DAPT:
    input:
        'output/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}.RData'
    output:
        'plot/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}_plots/FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_dlgs_plots_cont_DAPT_{pc}_{go}.txt'
    log:
        'logs/gsea_dlgs_plots_cont_DAPT_{pc}_{go}.log'
    shell:
        'src/gsea_dlgs_plots.sh {input} {output} >& {log}'

rule gsea_dlgs_plots_cont_DAPT_cov:
    input:
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_{go}.RData'
    output:
        'plot/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_{go}_plots/FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_dlgs_plots_cont_DAPT_cov_{pc}_{go}.txt'
    log:
        'logs/gsea_dlgs_plots_cont_DAPT_cov_{pc}_{go}.log'
    shell:
        'src/gsea_dlgs_plots.sh {input} {output} >& {log}'

#######################################
# DLGs Waterfall Plots
#######################################
rule dlgs_waterfall_cont_DAPT:
    input:
        'output/cont_DAPT/dlgs/{pc}/dlgs.tsv'
    output:
        'plot/cont_DAPT/dlgs/{pc}/dlgs_waterfall.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_waterfall_cont_DAPT_{pc}.txt'
    log:
        'logs/dlgs_waterfall_cont_DAPT_{pc}.log'
    shell:
        'src/dlgs_waterfall.sh {input} {output} >& {log}'

rule dlgs_waterfall_cont_DAPT_cov:
    input:
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs.tsv'
    output:
        'plot/cont_DAPT_cov/dlgs/{pc}/dlgs_waterfall.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_waterfall_cont_DAPT_cov_{pc}.txt'
    log:
        'logs/dlgs_waterfall_cont_DAPT_cov_{pc}.log'
    shell:
        'src/dlgs_waterfall.sh {input} {output} >& {log}'

#######################################
# DLGs GSEA Dot Plots
#######################################
rule dlgs_gsea_dotplot_cont_DAPT:
    input:
        'output/cont_DAPT/dlgs/{pc}/dlgs_gsea_bp.RData',
        'output/cont_DAPT/dlgs/{pc}/dlgs_gsea_mf.RData',
        'output/cont_DAPT/dlgs/{pc}/dlgs_gsea_cc.RData'
    output:
        'plot/cont_DAPT/dlgs/{pc}/dlgs_gsea_dotplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_gsea_dotplot_cont_DAPT_{pc}.txt'
    log:
        'logs/dlgs_gsea_dotplot_cont_DAPT_{pc}.log'
    shell:
        'src/dlgs_gsea_dotplot.sh {input} {output} >& {log}'

rule dlgs_gsea_dotplot_cont_DAPT_cov:
    input:
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_bp.RData',
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_mf.RData',
        'output/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_cc.RData'
    output:
        'plot/cont_DAPT_cov/dlgs/{pc}/dlgs_gsea_dotplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_gsea_dotplot_cont_DAPT_cov_{pc}.txt'
    log:
        'logs/dlgs_gsea_dotplot_cont_DAPT_cov_{pc}.log'
    shell:
        'src/dlgs_gsea_dotplot.sh {input} {output} >& {log}'

#######################################
# GSEA Basin (Oulhen)
#######################################
rule gsea_basin_Oulhen:
    input:
        'output/oulhen/{sample}/{pc}/seurat_landscaper.RData',
        'data/go_{go}_hpbase.RData'
    output:
        'plot/oulhen/{sample}/{pc}/basin_gsea_{go}.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in OULHEN_SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_basin_Oulhen_{sample}_{pc}_{go}.txt'
    log:
        'logs/gsea_basin_Oulhen_{sample}_{pc}_{go}.log'
    shell:
        'src/gsea_basin.sh {input} {output} >& {log}'

#######################################
# GSEA Basin Plots (Oulhen)
#######################################
rule gsea_basin_plots_Oulhen:
    input:
        'plot/oulhen/{sample}/{pc}/basin_gsea_{go}.RData'
    output:
        'plot/oulhen/{sample}/{pc}/basin_gsea_{go}_plots/FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in OULHEN_SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_basin_plots_Oulhen_{sample}_{pc}_{go}.txt'
    log:
        'logs/gsea_basin_plots_Oulhen_{sample}_{pc}_{go}.log'
    shell:
        'src/gsea_basin_plots.sh {input} {output} >& {log}'

#######################################
# GSEA for DLGs (Oulhen)
#######################################
rule gsea_dlgs_cont_DAPT_Oulhen:
    input:
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs.tsv',
        'data/go_{go}_hpbase.RData'
    output:
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}.RData'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_dlgs_cont_DAPT_Oulhen_{pc}_{go}.txt'
    log:
        'logs/gsea_dlgs_cont_DAPT_Oulhen_{pc}_{go}.log'
    shell:
        'src/gsea_dlgs.sh {input} {output} >& {log}'

#######################################
# GSEA Plots for DLGs (Oulhen)
#######################################
rule gsea_dlgs_plots_cont_DAPT_Oulhen:
    input:
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}.RData'
    output:
        'plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_{go}_plots/FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_dlgs_plots_cont_DAPT_Oulhen_{pc}_{go}.txt'
    log:
        'logs/gsea_dlgs_plots_cont_DAPT_Oulhen_{pc}_{go}.log'
    shell:
        'src/gsea_dlgs_plots.sh {input} {output} >& {log}'

#######################################
# DLGs Waterfall Plots (Oulhen)
#######################################
rule dlgs_waterfall_cont_DAPT_Oulhen:
    input:
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs.tsv'
    output:
        'plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_waterfall.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_waterfall_cont_DAPT_Oulhen_{pc}.txt'
    log:
        'logs/dlgs_waterfall_cont_DAPT_Oulhen_{pc}.log'
    shell:
        'src/dlgs_waterfall.sh {input} {output} >& {log}'

#######################################
# DLGs GSEA Dot Plots (Oulhen)
#######################################
rule dlgs_gsea_dotplot_cont_DAPT_Oulhen:
    input:
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_bp.RData',
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_mf.RData',
        'output/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_cc.RData'
    output:
        'plot/oulhen/cont_DAPT/dlgs/{pc}/dlgs_gsea_dotplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dlgs_gsea_dotplot_cont_DAPT_Oulhen_{pc}.txt'
    log:
        'logs/dlgs_gsea_dotplot_cont_DAPT_Oulhen_{pc}.log'
    shell:
        'src/dlgs_gsea_dotplot.sh {input} {output} >& {log}'

#######################################
# GSEA Basin Hexbin Plots
#######################################
rule gsea_basin_hexbin_plots:
    input:
        'plot/{sample}/{pc}/basin_gsea_{go}.RData',
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/basin_gsea_{go}_hexbin/FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_basin_hexbin_plots_{sample}_{pc}_{go}.txt'
    log:
        'logs/gsea_basin_hexbin_plots_{sample}_{pc}_{go}.log'
    shell:
        'src/gsea_basin_hexbin_plots.sh {input} {output} >& {log}'

#######################################
# GSEA Basin Hexbin Plots (Oulhen)
#######################################
rule gsea_basin_hexbin_plots_Oulhen:
    input:
        'plot/oulhen/{sample}/{pc}/basin_gsea_{go}.RData',
        'output/oulhen/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/oulhen/{sample}/{pc}/basin_gsea_{go}_hexbin/FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in OULHEN_SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/gsea_basin_hexbin_plots_Oulhen_{sample}_{pc}_{go}.txt'
    log:
        'logs/gsea_basin_hexbin_plots_Oulhen_{sample}_{pc}_{go}.log'
    shell:
        'src/gsea_basin_hexbin_plots.sh {input} {output} >& {log}'

