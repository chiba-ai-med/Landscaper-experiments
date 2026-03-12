import pandas as pd
import re
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

PCS = list(range(3, 13))
SAMPLES = ['cont', 'DAPT', 'integrated', 'cont_cov', 'DAPT_cov', 'integrated_cov']
CELLTYPES = ['Aboral_ectoderm', 'Oral_ectoderm', 'Ciliary_band', 'Anus', 'Neurons']
BASIN_STATUS = ['basin', 'nonbasin']

container: 'docker://koki/urchin_workflow_seurat:20251014'

rule all:
    input:
        expand('plot/{sample}/{pc}/P_metropolis.tsv', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber.tsv', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_metropolis_rw.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_glauber_rw.RData', sample=SAMPLES, pc=PCS),
        expand('plot/{sample}/{pc}/P_metropolis_rw_{ct}_{bs}.png',
            sample=SAMPLES, pc=PCS, ct=CELLTYPES, bs=BASIN_STATUS),
        expand('plot/{sample}/{pc}/P_glauber_rw_{ct}_{bs}.png',
            sample=SAMPLES, pc=PCS, ct=CELLTYPES, bs=BASIN_STATUS)

#######################################
# Transition Probability
#######################################
rule transition_probability:
    input:
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/P_metropolis.tsv',
        'plot/{sample}/{pc}/P_glauber.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/transition_probability_{sample}_{pc}.txt'
    log:
        'logs/transition_probability_{sample}_{pc}.log'
    shell:
        'src/transition_probability.sh {input} {output} >& {log}'

#######################################
# Random Walk
#######################################
rule random_walk:
    input:
        'plot/{sample}/{pc}/P_metropolis.tsv',
        'plot/{sample}/{pc}/P_glauber.tsv'
    output:
        'plot/{sample}/{pc}/P_metropolis_rw.RData',
        'plot/{sample}/{pc}/P_glauber_rw.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/random_walk_{sample}_{pc}.txt'
    log:
        'logs/random_walk_{sample}_{pc}.log'
    shell:
        'src/random_walk.sh {input} {output} >& {log}'

#######################################
# Plot Random Walk
#######################################
rule plot_random_walk_metropolis:
    input:
        'output/integrated/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_metropolis_rw.RData',
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/P_metropolis_rw_Aboral_ectoderm_basin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Oral_ectoderm_basin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Ciliary_band_basin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Anus_basin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Neurons_basin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Aboral_ectoderm_nonbasin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Oral_ectoderm_nonbasin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Ciliary_band_nonbasin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Anus_nonbasin.png',
        'plot/{sample}/{pc}/P_metropolis_rw_Neurons_nonbasin.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_metropolis_{sample}_{pc}.txt'
    log:
        'logs/plot_random_walk_metropolis_{sample}_{pc}.log'
    shell:
        'src/plot_random_walk.sh {input} {output} >& {log}'

rule plot_random_walk_glauber:
    input:
        'output/integrated/{pc}/seurat_landscaper.RData',
        'plot/{sample}/{pc}/P_glauber_rw.RData',
        'output/{sample}/{pc}/seurat_landscaper.RData'
    output:
        'plot/{sample}/{pc}/P_glauber_rw_Aboral_ectoderm_basin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Oral_ectoderm_basin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Ciliary_band_basin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Anus_basin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Neurons_basin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Aboral_ectoderm_nonbasin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Oral_ectoderm_nonbasin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Ciliary_band_nonbasin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Anus_nonbasin.png',
        'plot/{sample}/{pc}/P_glauber_rw_Neurons_nonbasin.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_glauber_{sample}_{pc}.txt'
    log:
        'logs/plot_random_walk_glauber_{sample}_{pc}.log'
    shell:
        'src/plot_random_walk.sh {input} {output} >& {log}'
