# Landscaper-experiments

Energy Landscape Analysis workflows for single-cell transcriptomics using [Landscaper](https://github.com/chiba-ai-med/Landscaper).

This workflow consists of 10 workflows as follows:

- **workflow/preprocess.smk**: PCA, binarization, and sample stratification

- **workflow/landscaper.smk**: Landscaper execution, postprocessing, and DLGs computation

- **workflow/rank_estimate.smk**: PC rank selection by AUC cross-validation

- **workflow/random_walk.smk**: Transition probability matrices (Metropolis/Glauber) and random walk visualization

- **workflow/vector_field.smk**: Coarse-grained vector fields

- **workflow/graph_embedding.smk**: Graph embedding of state transition networks

- **workflow/discon_graph.smk**: Disconnectivity graphs and dendrograms

- **workflow/fate_probability.smk**: Fate probability (absorption probabilities)

- **workflow/gsea.smk**: GSEA for Basin and DLGs (waterfall plot, dot plot)

- **workflow/oulhen.smk**: Oulhen dataset pipeline (encompasses most of the above, plus cross-species DLGs overlap analysis)

## Requirements

- Bash: GNU bash, version 4.2.46(1)-release (x86_64-redhat-linux-gnu)
- Snakemake: 6.5.3
- Singularity: 3.8.0

## How to reproduce this workflow

### In Local Machine

```
snakemake -s workflow/preprocess.smk -j 4 --use-singularity
snakemake -s workflow/landscaper.smk -j 4 --use-singularity
snakemake -s workflow/rank_estimate.smk -j 4 --use-singularity
snakemake -s workflow/random_walk.smk -j 4 --use-singularity
snakemake -s workflow/vector_field.smk -j 4 --use-singularity
snakemake -s workflow/graph_embedding.smk -j 4 --use-singularity
snakemake -s workflow/discon_graph.smk -j 4 --use-singularity
snakemake -s workflow/fate_probability.smk -j 4 --use-singularity
snakemake -s workflow/gsea.smk -j 4 --use-singularity
snakemake -s workflow/oulhen.smk -j 4 --use-singularity
```

### In Open Grid Engine

```
snakemake -s workflow/preprocess.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/landscaper.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/rank_estimate.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/random_walk.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/vector_field.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/graph_embedding.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/discon_graph.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/fate_probability.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/gsea.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/oulhen.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
```

### In Slurm

```
snakemake -s workflow/preprocess.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/landscaper.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/rank_estimate.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/random_walk.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/vector_field.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/graph_embedding.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/discon_graph.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/fate_probability.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/gsea.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/oulhen.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
```

## Authors

- Koki Tsuyuzaki
