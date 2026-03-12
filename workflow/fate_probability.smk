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
        expand('plot/{sample}/{pc}/P_metropolis_fp.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_fp.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_metropolis_fp/FINISH', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_fp/FINISH', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_metropolis_fp/dimplot_FINISH', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_fp/dimplot_FINISH', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_metropolis_fp/scatter_FINISH', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_fp/scatter_FINISH', sample=SAMPLES, pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont.RData', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont.RData', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT.RData', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT.RData', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont/FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont/FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT/FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT/FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont/dimplot_FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont/dimplot_FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT/dimplot_FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT/dimplot_FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont/scatter_FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont/scatter_FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT/scatter_FINISH', pc=PCS),
        expand('plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT/scatter_FINISH', pc=PCS)

#######################################
# Fate Probability
#######################################
rule fate_probability:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis.tsv',
        'plot/{sample}/{pc}/P_glauber.tsv'
    output:
        'plot/{sample}/{pc}/P_metropolis_fp.RData',
        'plot/{sample}/{pc}/P_glauber_fp.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/fate_probability_{sample}_{pc}.txt'
    log:
        'logs/fate_probability_{sample}_{pc}.log'
    shell:
        'src/fate_probability.sh {input} {output} >& {log}'

#######################################
# Heatmap Fate Probability
#######################################
rule heatmap_fate_probability:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis_fp.RData',
        'plot/{sample}/{pc}/P_glauber_fp.RData'
    output:
        'plot/{sample}/{pc}/P_metropolis_fp/FINISH',
        'plot/{sample}/{pc}/P_glauber_fp/FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/heatmap_fate_probability_{sample}_{pc}.txt'
    log:
        'logs/heatmap_fate_probability_{sample}_{pc}.log'
    shell:
        'src/heatmap_fate_probability.sh {input} {output} >& {log}'

#######################################
# Dimplot Fate Probability
#######################################
rule dimplot_fate_probability:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis_fp.RData',
        'plot/{sample}/{pc}/P_glauber_fp.RData'
    output:
        'plot/{sample}/{pc}/P_metropolis_fp/dimplot_FINISH',
        'plot/{sample}/{pc}/P_glauber_fp/dimplot_FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_fate_probability_{sample}_{pc}.txt'
    log:
        'logs/dimplot_fate_probability_{sample}_{pc}.log'
    shell:
        'src/dimplot_fate_probability.sh {input} {output} >& {log}'

#######################################
# Energy Y-axis Range (per PC)
#######################################
rule energy_ylim:
    input:
        expand('plot/{sample}/{{pc}}/Landscaper/E.tsv', sample=SAMPLES)
    output:
        'output/{pc}/energy_ylim.tsv'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/energy_ylim_{pc}.txt'
    log:
        'logs/energy_ylim_{pc}.log'
    shell:
        'src/energy_ylim.sh {input} {output} >& {log}'

#######################################
# Scatter Fate Probability vs Energy
#######################################
rule scatter_fate_energy_metropolis:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis_fp.RData',
        'plot/{sample}/{pc}/Landscaper/igraph.RData',
        'output/{pc}/energy_ylim.tsv'
    output:
        'plot/{sample}/{pc}/P_metropolis_fp/scatter_FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/scatter_fate_energy_metropolis_{sample}_{pc}.txt'
    log:
        'logs/scatter_fate_energy_metropolis_{sample}_{pc}.log'
    shell:
        'src/scatter_fate_energy.sh {input} {output} >& {log}'

rule scatter_fate_energy_glauber:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_glauber_fp.RData',
        'plot/{sample}/{pc}/Landscaper/igraph.RData',
        'output/{pc}/energy_ylim.tsv'
    output:
        'plot/{sample}/{pc}/P_glauber_fp/scatter_FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/scatter_fate_energy_glauber_{sample}_{pc}.txt'
    log:
        'logs/scatter_fate_energy_glauber_{sample}_{pc}.log'
    shell:
        'src/scatter_fate_energy.sh {input} {output} >& {log}'

#######################################
# Fate Probability (cont_DAPT_cov)
#######################################
rule fate_probability_cont_DAPT_cov:
    input:
        'output/cont_cov/{pc}/seurat_landscaper.RData',
        'output/DAPT_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_cov/{pc}/P_metropolis.tsv',
        'plot/DAPT_cov/{pc}/P_metropolis.tsv',
        'plot/cont_cov/{pc}/P_glauber.tsv',
        'plot/DAPT_cov/{pc}/P_glauber.tsv'
    output:
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont.RData',
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT.RData'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/fate_probability_cont_DAPT_cov_{pc}.txt'
    log:
        'logs/fate_probability_cont_DAPT_cov_{pc}.log'
    shell:
        'src/fate_probability_cont_DAPT_cov.sh {input} {output} >& {log}'

#######################################
# Heatmap Fate Probability (cont_cov)
#######################################
rule heatmap_fate_probability_cont_cov:
    input:
        'output/cont_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont.RData'
    output:
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont/FINISH',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont/FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/heatmap_fate_probability_cont_cov_{pc}.txt'
    log:
        'logs/heatmap_fate_probability_cont_cov_{pc}.log'
    shell:
        'src/heatmap_fate_probability.sh {input} {output} >& {log}'

#######################################
# Heatmap Fate Probability (DAPT_cov)
#######################################
rule heatmap_fate_probability_DAPT_cov:
    input:
        'output/DAPT_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT.RData'
    output:
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT/FINISH',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT/FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/heatmap_fate_probability_DAPT_cov_{pc}.txt'
    log:
        'logs/heatmap_fate_probability_DAPT_cov_{pc}.log'
    shell:
        'src/heatmap_fate_probability.sh {input} {output} >& {log}'

#######################################
# Dimplot Fate Probability (cont_cov)
#######################################
rule dimplot_fate_probability_cont_cov:
    input:
        'output/cont_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont.RData'
    output:
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont/dimplot_FINISH',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont/dimplot_FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_fate_probability_cont_cov_{pc}.txt'
    log:
        'logs/dimplot_fate_probability_cont_cov_{pc}.log'
    shell:
        'src/dimplot_fate_probability.sh {input} {output} >& {log}'

#######################################
# Dimplot Fate Probability (DAPT_cov)
#######################################
rule dimplot_fate_probability_DAPT_cov:
    input:
        'output/DAPT_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT.RData'
    output:
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT/dimplot_FINISH',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT/dimplot_FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_fate_probability_DAPT_cov_{pc}.txt'
    log:
        'logs/dimplot_fate_probability_DAPT_cov_{pc}.log'
    shell:
        'src/dimplot_fate_probability.sh {input} {output} >& {log}'

#######################################
# Scatter Fate vs Energy (cont_cov)
#######################################
rule scatter_fate_energy_cont_cov_metropolis:
    input:
        'output/cont_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont.RData',
        'plot/cont_cov/{pc}/Landscaper/igraph.RData',
        'output/{pc}/energy_ylim.tsv'
    output:
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_cont/scatter_FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/scatter_fate_energy_cont_cov_metropolis_{pc}.txt'
    log:
        'logs/scatter_fate_energy_cont_cov_metropolis_{pc}.log'
    shell:
        'src/scatter_fate_energy.sh {input} {output} >& {log}'

rule scatter_fate_energy_cont_cov_glauber:
    input:
        'output/cont_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont.RData',
        'plot/cont_cov/{pc}/Landscaper/igraph.RData',
        'output/{pc}/energy_ylim.tsv'
    output:
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_cont/scatter_FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/scatter_fate_energy_cont_cov_glauber_{pc}.txt'
    log:
        'logs/scatter_fate_energy_cont_cov_glauber_{pc}.log'
    shell:
        'src/scatter_fate_energy.sh {input} {output} >& {log}'

#######################################
# Scatter Fate vs Energy (DAPT_cov)
#######################################
rule scatter_fate_energy_DAPT_cov_metropolis:
    input:
        'output/DAPT_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT.RData',
        'plot/DAPT_cov/{pc}/Landscaper/igraph.RData',
        'output/{pc}/energy_ylim.tsv'
    output:
        'plot/cont_DAPT_cov/{pc}/P_metropolis_fp_DAPT/scatter_FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/scatter_fate_energy_DAPT_cov_metropolis_{pc}.txt'
    log:
        'logs/scatter_fate_energy_DAPT_cov_metropolis_{pc}.log'
    shell:
        'src/scatter_fate_energy.sh {input} {output} >& {log}'

rule scatter_fate_energy_DAPT_cov_glauber:
    input:
        'output/DAPT_cov/{pc}/seurat_landscaper.RData',
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT.RData',
        'plot/DAPT_cov/{pc}/Landscaper/igraph.RData',
        'output/{pc}/energy_ylim.tsv'
    output:
        'plot/cont_DAPT_cov/{pc}/P_glauber_fp_DAPT/scatter_FINISH'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/scatter_fate_energy_DAPT_cov_glauber_{pc}.txt'
    log:
        'logs/scatter_fate_energy_DAPT_cov_glauber_{pc}.log'
    shell:
        'src/scatter_fate_energy.sh {input} {output} >& {log}'
